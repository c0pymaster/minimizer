// #define NDEBUG // uncomment to turn off assert checks (increases speed)

#include "interval.hpp"
#include "geometry.h"
#include "melzak.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <utility>
#include <assert.h>
#include <future>
#include <chrono>
#include <queue>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <set>

#define eprintf(...) fprintf(stderr, __VA_ARGS__), fflush(stderr)
#define sz(a) ((int) (a).size())

using namespace std;
using namespace cxsc;

// usage: ./main threads
// Parameter threads is the number of threads used. If threads = 0, the program tries
// to guess the number of concurrent threads supported by the machine.
//
// The program verifies that L(p) >= L_0 for every point p = (x, y, alpha, xi, xi_1, xi_2) in P \ P_0.
//
// It is assumed that r = 1 (all the points can be scaled by 1 / r, if it is not so).
//
// Tasks are run in parallel using multiple threads.
// Each task corresponds to a box, the initial task corresponds to the whole parametric space.
//
// For each task the program either verifies the inequality for all points in the corresponding box,
// or splits the box into two smaller boxes and considers each of them as a separate task.
// The box (task) of depth i is the box which was obtained from the initial task by i splits.
//
// Interval arithmetics are used to ensure that all computations are correct
// (see http://www.xsc.de for information on the library used).

const int maxDepth = 90;

// constants used in error estimation
const interval C1 = 2 * (cos(11 * pi / 24) + (cos(11 * pi / 8) + sin(11 * pi / 24)) / sin60);
const interval C2 = power(3 * cos(17 * pi / 24), 2) + power(cos(13 * pi / 8), 2);
const interval C3 = 6 * cos(17 * pi / 24) * cos(13 * pi / 8);

// dimensions of P_0
const interval radXY = interval(1) / interval(10), radAlpha = interval(1) / interval(20);

// sum of alpha and beta
const interval sumAlphaBeta = pi * 11 / 12;

// d = |A_1C_1| = |A_1C_2|
const interval d = interval(4) / sqrt(interval(6)) + interval(4);

// part of the cycle: Z_1-W_1-V-W_2-Z_2
struct cyclePart {
    // lines through segments Z_1W_1 and Z_2W_2;
    // Z_1W_1.o = W_1, Z_2W_2.o = W_2
    line Z_1W_1, Z_2W_2;
    point V;

    // total length of Z_1-W_1-V-W_2-Z_2
    interval len;

    cyclePart() {}

    // given x, y, alpha, find lines Z_1W_1, Z_2W_2 and point V
    cyclePart(const interval & x, const interval & y, const interval & alpha) {
        interval beta = sumAlphaBeta - alpha;

        // W_1 and W_2
        Z_1W_1.o = point(x + cos(alpha), sin(alpha));
        Z_2W_2.o = point(sin(beta), y + cos(beta));
        
        Z_1W_1.v = point(cos(2 * alpha), sin(2 * alpha));
        Z_2W_2.v = point(sin(2 * beta), cos(2 * beta));
        V = intersect(Z_1W_1, Z_2W_2);

        len = 2 * d - Z_1W_1.o.x - Z_2W_2.o.y + (V - Z_1W_1.o).len() + (V - Z_2W_2.o).len();
    }
};

// x_0, y_0, alpha_0, and L_0
const interval x_0 = sqrt(interval(2)), y_0 = sqrt(interval(2)), alpha_0 = pi * 11 / 24;
const cyclePart C_0(x_0, x_0, alpha_0);
const interval L_0 = C_0.len + C_0.V.len() - 1;

// time point of the start of the program
chrono::steady_clock::time_point start = chrono::steady_clock::now();

// return the amount of seconds elapsed since the start of the program
int getSecondsElapsed() {
    return (int) chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - start).count();
}

// format seconds to hh:mm:ss and print to the terminal
void eprintTime(int seconds) {
    int minutes = seconds / 60, hours = minutes / 60;
    seconds %= 60, minutes %= 60;
    eprintf("%4d:%02d:%02d", hours, minutes, seconds);
}

// a point in the parametric space
struct point6d {
    interval x, y, alpha, xi, xi_1, xi_2;

    point6d(): x(0), y(0), alpha(0), xi(0), xi_1(0), xi_2(0) {}
    point6d(const interval & x, const interval & y, const interval & alpha,
            const interval & xi, const interval & xi_1, const interval & xi_2):
                x(x), y(y), alpha(alpha), xi(xi), xi_1(xi_1), xi_2(xi_2) {}

    // return the point, multiplied by the constant
    point6d operator *(const interval & m) const {
        return point6d(x * m, y * m, alpha * m, xi * m, xi_1 * m, xi_2 * m);
    }

    // return the sum of two points
    point6d operator +(const point6d & p) const {
        return point6d(x + p.x, y + p.y, alpha + p.alpha, xi + p.xi, xi_1 + p.xi_1, xi_2 + p.xi_2);
    }
};

ostream & operator <<(ostream & os, const point6d & p) {
   os << p.x << ' ' << p.y << ' ' << p.alpha << ' ' << p.xi << ' ' << p.xi_1 << ' ' << p.xi_2;
   return os;
}

istream & operator >>(istream & is, point6d & p) {
    is >> p.x >> p.y >> p.alpha >> p.xi >> p.xi_1 >> p.xi_2;
    return is;
}

// a task corresponding to a box P(c, delta)
struct Task {
    point6d c, delta;
    int depth;

    Task(): c(), delta(), depth() {}
    Task(const point6d & c, const point6d & delta, const int & depth): c(c), delta(delta), depth(depth) {}

    bool operator <(const Task & t) const {
        return depth > t.depth;
    }
};

struct Verifier {
    int threads;

    // thread pool maintains a priority queue of tasks and distributes them among threads;
    // tasks of higher depth have higher priority
    struct ThreadPool {
        vector<thread> threadPool;

        // priority queue of tasks
        multiset<Task> tasks;

        // tasksProcessed[i] is the amount of processed tasks of depth i
        atomic_llong tasksProcessed[maxDepth];
        // tasksToProcess[i] is the amount of tasks added to the queue but not yet processed
        atomic_int tasksToProcess[maxDepth];

        // lastTaskTaken[i] is the time point of the latest moment when a task of depth i has been taken
        atomic_int lastTaskTaken[maxDepth];

        mutex tasksMutex;
        condition_variable condition;
        atomic_bool terminate;
        atomic_bool stopped;
        
        // initialize threads
        void init(int threads, Verifier *verifier) {
            terminate = false;
            stopped = false;

            for (int i = 0; i < threads; i++) {
                threadPool.emplace_back(thread(&ThreadPool::Runner, this, verifier));
            }
        }

        ~ThreadPool() {
            if (!stopped.load()) {
                ShutDown();
            }
        }

        // add a new task t
        void addTask(const Task & t) {
            {
                unique_lock<mutex> lock(tasksMutex);
                tasks.insert(t);
                ++tasksToProcess[t.depth];
            }

            // wake up one (random) thread.
            condition.notify_one();
        }

        // terminate all threads
        void ShutDown() {
            terminate = true;
            condition.notify_all();

            for (thread &thread : threadPool) {
                thread.join();
            }

            threadPool.clear();
            stopped = true;
        }

        // the function run by each thread; takes a task of highest available depth and runs it;
        // idle if no tasks are pending
        void Runner(Verifier *verifier) {
            Task t;
            Melzak melzak;

            while (true) {
                {
                    unique_lock<mutex> lock(tasksMutex);

                    condition.wait(lock, [this] { return !tasks.empty() || terminate.load(); });

                    if (terminate.load()) {
                        return;
                    }

                    t = *tasks.begin();
                    tasks.erase(tasks.begin());
                }

                lastTaskTaken[t.depth] = getSecondsElapsed();

                // run the chosen task
                verifier->checkBox(t, melzak);
                --tasksToProcess[t.depth];
                ++tasksProcessed[t.depth];
            }
        }
    } pool;

    mutex errorMutex;
    condition_variable errorCondition;
    atomic_bool errorOccured;
    string errorMessage;

    // once a second show current progress in the terminal
    void show() {
        for (int wakeUp = 1; ; ++wakeUp) {
            {
                unique_lock<mutex> lock(errorMutex);
                // wait until a second passes or an error occurs
                errorCondition.wait_until(lock, start + std::chrono::seconds(wakeUp), [this] { return errorOccured.load(); });
            }

            // time elapsed since the start of the program
            int now = getSecondsElapsed();

            eprintf("\33[2J\33[H");
            eprintf("time elapsed: ");
            eprintTime(now + saveTime);
            eprintf("\n");

            // stop if an error occured; show error message in the terminal
            if (errorOccured.load()) {
                pool.ShutDown();
                eprintf("Error occured: %s\n", errorMessage.c_str());
                return;
            }

            bool noTasksLeft = true;

            for (int _ = 0; _ < 3; ++_) {
                eprintf(" |       |         task(s) | task(s) |      time since");
            }
            eprintf(" |\n");
            for (int _ = 0; _ < 3; ++_) {
                eprintf(" | depth |       processed | pending | last task taken");
            };
            eprintf(" |\n");
            for (int _ = 0; _ < 3; ++_) {
                eprintf(" | ----- | --------------- | ------- | ---------------");
            };

            for (int it = 0; it < maxDepth; ++it) {
                int i = maxDepth / 3 * (it % 3) + it / 3;
                if (pool.tasksProcessed[i].load() != 0 || pool.tasksToProcess[i].load() != 0) {
                    if (it % 3 == 0) {
                        eprintf(" |\n");
                    }
                    eprintf(" | %5d | %15lld | %7d |      ", i, pool.tasksProcessed[i].load(), pool.tasksToProcess[i].load());
                    if (pool.tasksToProcess[i].load() != 0) {
                        noTasksLeft = false;
                        eprintTime(now - pool.lastTaskTaken[i].load());
                    } else {
                        eprintf("          ");
                    }
                }
            }

            // once a minute save current progress to the file
            if (wakeUp % 60 == 0) {
                eprintf("\nSaving...\n");
                save(now);
            }

            // stop if there are no tasks left
            if (noTasksLeft) {
                save(now);
                eprintf("\nFinished!\n");
                break;
            }
        }
    }

    // time spent in previous runs of the program (zero if there is no save file)
    int saveTime = 0;

    // save current progress to the file
    void save(int now) {
        pool.ShutDown();

        ofstream saveFile;
        saveFile.open("save.txt", ios_base::out);
        saveFile << SetPrecision(18, 15);
        saveFile << now + saveTime << '\n';
        for (int i = 0; i < maxDepth; ++i) {
            saveFile << pool.tasksProcessed[i].load() << ' ';
        }
        saveFile << '\n';
        for (int i = 0; i < maxDepth; ++i) {
            saveFile << pool.lastTaskTaken[i].load() - now << ' ';
        }
        saveFile << '\n';
        for (auto t : pool.tasks) {
            saveFile << t.c << ' ' << t.delta << ' ' << t.depth << '\n';
        }
        saveFile.close();

        pool.init(threads, this);
    }

    // initialize verifier
    void init(char** argv) {
        sscanf(argv[1], "%d", &threads);
        if (threads == 0) {
            // try to guess the number of concurrent threads supported by the machine
            threads = thread::hardware_concurrency() - 1;
        }
        
        errorOccured = false;

        // initialize thread pool
        pool.init(threads, this);
    }

    // consider the box P(t.c, t.delta)
    // if L(c) - Err(c, delta) < L_0, divide the box into 2 smaller boxes
    // by splitting it in the dimension corresponding to the largest error term,
    // and add a task for each of them
    void checkBox(const Task & t, Melzak & melzak) {
        const point6d & c = t.c;
        const point6d & delta = t.delta;
        
        if (getSign(radXY - delta.x - abs(c.x - x_0)) == 1 &&
            getSign(radXY - delta.y - abs(c.y - y_0)) == 1 &&
            getSign(radAlpha - delta.alpha - abs(c.alpha - alpha_0)) == 1) {
            // P(c, delta) is contained in P_0
            return;
        }
        if (getSign(4 - power(c.x + delta.x, 2) - power(c.y + delta.y, 2)) == 1) {
            // B_r(y(W_1)) intersects B_r(y(W_2)) for each point in the box
            return;
        }
      
        cyclePart cycle(c.x, c.y, c.alpha);

        interval beta = sumAlphaBeta - c.alpha;

        point p[4];
        p[0] = cycle.V;
        p[1] = point(cos(c.xi), sin(c.xi));
        p[2] = point(c.x - cos(c.xi_1), sin(c.xi_1));
        p[3] = point(sin(c.xi_2), c.y - cos(c.xi_2));

        // length = L(c)
        interval length = melzak.solve(p) + cycle.len;

        if (getSign(L_0 - length) == 1) {
            // construction shorter than Sigma_0 found; stop with an error
            {
                unique_lock<mutex> lock(errorMutex);
                ostringstream ss;
                ss << "shorter construction found: " << length << " < " << L_0 << '\n';
                ss << "Steiner part length: " << length - cycle.len << "\n";
                ss << " V: " << p[0] << '\n';
                ss << " Q: " << p[1] << '\n';
                ss << "Q1: " << p[2] << '\n';
                ss << "Q2: " << p[3] << '\n';
                ss << "configuration point: " << c << '\n';
                errorMessage = ss.str();
                errorOccured = true;
            }

            errorCondition.notify_all();
            return;
        }

        // calculate Err(c, delta)
        point6d err;

        err.alpha = 2 * sqrt(power(c.x + delta.x, 2) + power(c.y + delta.y, 2));
        err.alpha += sqrt(C2 - C3 * cos(2 * (abs(c.alpha - 11 * pi / 24) + delta.alpha)));
        err.alpha /= sin60;
        
        err.x = sin(2 * (c.alpha - delta.alpha)) / sin60;
        err.y = sin(2 * (beta - delta.alpha)) / sin60;

        err.x += 1, err.y += 1;
        err.xi = 1, err.xi_1 = 1, err.xi_2 = 1;

        interval tmpSinA = sin(2 * (c.alpha - delta.alpha - pi / 3)), tmpSinB = sin(2 * (c.alpha + delta.alpha - pi / 3));
        interval tmpCosA = cos(2 * (c.alpha - delta.alpha - pi / 3)), tmpCosB = cos(2 * (c.alpha + delta.alpha - pi / 3));

        err.x += 1 - tmpSinA / sin60;
        err.y += 1 - tmpCosB / sin60;

        interval F1 =  (c.x + delta.x) * tmpCosA - (c.y - delta.y) * tmpSinA;
        F1 = F1 * 2 / sin60 + C1 * sin((c.alpha - beta) / 2 + delta.alpha);
        interval F2 = -(c.x - delta.x) * tmpCosB + (c.y + delta.y) * tmpSinB;
        F2 = F2 * 2 / sin60 - C1 * sin((c.alpha - beta) / 2 - delta.alpha);
        
        err.alpha += mymax(interval(0), mymax(F1, F2));

        err.x *= delta.x, err.y *= delta.y, err.alpha *= delta.alpha;
        err.xi *= delta.xi, err.xi_1 *= delta.xi_1, err.xi_2 *= delta.xi_2;

        interval Err = err.x + err.y + err.alpha + err.xi + err.xi_1 + err.xi_2;
        interval mx = mymax(mymax(err.x, mymax(err.y, err.alpha)), mymax(err.xi, mymax(err.xi_1, err.xi_2)));
        
        if (getSign(length - Err - L_0) == 1) {
            // the inequality is verified for every point in the box
            return;
        }

        if (t.depth + 1 == maxDepth) {
            // maximum allowed depth exceeded; stop with an error
            {
                unique_lock<mutex> lock(errorMutex);
                errorMessage = "maximum depth exceeded\n";
                errorOccured = true;
            }

            errorCondition.notify_all();
            return;
        }

        // add new tasks
        point6d ndelta = delta, shift;
        // choose the dimension corresponding to the largest error term
        if (getSign(mx - err.x) == 0) {
            ndelta.x /= 2;
            shift.x = ndelta.x;
        } else if (getSign(mx - err.y) == 0) {
            ndelta.y /= 2;
            shift.y = ndelta.y;
        } else if (getSign(mx - err.alpha) == 0) {
            ndelta.alpha /= 2;
            shift.alpha = ndelta.alpha;
        } else if (getSign(mx - err.xi) == 0) {
            ndelta.xi /= 2;
            shift.xi = ndelta.xi;
        } else if (getSign(mx - err.xi_1) == 0) {
            ndelta.xi_1 /= 2;
            shift.xi_1 = ndelta.xi_1;
        } else {
            assert(getSign(mx - err.xi_2) == 0);
            ndelta.xi_2 /= 2;
            shift.xi_2 = ndelta.xi_2;
        }

        pool.addTask(Task(c + shift, ndelta, t.depth + 1));
        pool.addTask(Task(c + shift * interval(-1), ndelta, t.depth + 1));
    }

    // if there is no save file detected, add the initial task;
    // if there is a save file, load the progress from it
    void addTasks() {
        ifstream saveFile("save.txt");
        if (saveFile.is_open()) {
            string line;
            {
                getline(saveFile, line);
                istringstream iss(line);
                iss >> saveTime;
            }
            {
                getline(saveFile, line);
                istringstream iss(line);
                for (int i = 0; i < maxDepth; ++i) {
                    long long tmp;
                    iss >> tmp;
                    pool.tasksProcessed[i] = tmp;
                }
            }
            {
                getline(saveFile, line);
                istringstream iss(line);
                for (int i = 0; i < maxDepth; ++i) {
                    long long tmp;
                    iss >> tmp;
                    pool.lastTaskTaken[i] = tmp;
                }
            }
            while (getline(saveFile, line)) {
                istringstream iss(line);
                Task t;
                iss >> t.c >> t.delta >> t.depth;
                pool.addTask(t);
            }
            saveFile.close();
        } else {
            for (int i = 0; i < maxDepth; ++i) {
                pool.lastTaskTaken[i] = 0;
            }

            Task t;
            t.delta = point6d(d - 1, d - 1, pi / 12, pi / 2, pi / 2, pi / 2) * cos60;
            t.c = t.delta;
            t.c.alpha += sumAlphaBeta - pi / 2;
            t.depth = 0;
            pool.addTask(t);
        }
    }
} verifier;


int main(int argc, char** argv) {
    if (argc != 2) {
        eprintf("Wrong number of arguments\n");
        exit(1);
    }

    verifier.init(argv);
    verifier.addTasks();
    verifier.show();
    return 0;
}
