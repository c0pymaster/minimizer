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

const int maxDepth = 20;

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
// or splits the box into 64 smaller boxes and considers each of them as a separate task.
// The box (task) of depth i is the box which was obtained from the initial task by i splits.
//
// Interval arithmetics are used to ensure that all computations are correct
// (see http://www.xsc.de for information on the library used).

// a point in the parametric space
struct point6d {
    interval x, y, alpha, xi, xi_1, xi_2;

    point6d() {}
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

struct Verifier {
    // dimensions of P_0
    interval radXY = interval(1) / interval(10), radAlpha = interval(1) / interval(30);
    
    // sum of alpha and beta
    interval sumAlphaBeta = pi * 11 / 12;

    // d = |A_1C_1| = |A_1C_2|
    interval d = interval(4) / sqrt(interval(6)) + interval(4);

    int threads;

    // x_0, y_0, alpha_0, and L_0 (L_0 is calculated in init method)
    interval x_0 = sqrt(interval(2)), y_0 = sqrt(interval(2)), alpha_0 = pi * 11 / 24, L_0;

    // deltas[i] stores the delta values for the box of depth i
    vector<point6d> deltas;

    // shift[i] is the list of the vectors from the center of the box of depth i
    // to the centers of 64 smaller boxes into which it is divided;
    vector<vector<point6d> > shift;

    // thread pool maintains a stack of tasks for each depth and distributes them among threads;
    // tasks of higher depth have higher priority
    struct ThreadPool {
        vector<thread> threadPool;

        // tasks[i] is the stack of tasks of depth i
        vector<vector<point6d> > tasks;

        // total amount of tasks in all stacks
        int tasksAmount;

        mutex tasksMutex;
        condition_variable condition;
        atomic_bool terminate;
        atomic_bool stopped;
        
        // initialize threads
        void init(int threads, Verifier *verifier) {
            terminate = false;
            stopped = false;

            tasksAmount = 0;

            for (int i = 0; i < threads; i++) {
                threadPool.emplace_back(thread(&ThreadPool::Runner, this, verifier));
            }
        }

        ~ThreadPool() {
            if (!stopped) {
                ShutDown();
            }
        }

        // add a new task of depth d with the center at t
        void addTask(const point6d & t, int d) {
            {
                unique_lock<mutex> lock(tasksMutex);
                for (; sz(tasks) <= d; tasks.push_back(vector<point6d>()));
                tasks[d].push_back(t);
                ++tasksAmount;
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

        // the function run by each thread; takes a task with highest available depth and runs it;
        // idle if no tasks are pending
        void Runner(Verifier *verifier) {
            point6d t;
            int depth;
            Melzak melzak;

            while (true) {
                {
                    unique_lock<mutex> lock(tasksMutex);

                    condition.wait(lock, [this] { return tasksAmount > 0 || terminate; });

                    if (terminate) {
                        return;
                    }

                    for (int i = sz(tasks) - 1; i >= 0; --i) {
                        if (sz(tasks[i]) > 0) {
                            depth = i;
                            t = *tasks[i].rbegin();
                            tasks[i].pop_back();
                            --tasksAmount;
                            break;
                        }
                    }
                }

                {
                    unique_lock<mutex> lock(verifier->timeMutex);
                    verifier->lastTaskTaken[depth] = chrono::steady_clock::now();
                }

                // run the chosen task
                verifier->checkBox(t, melzak, depth);
                --verifier->tasksToProcess[depth];
                ++verifier->tasksProcessed[depth];
            }
        }
    } pool;

    // tasksProcessed[i] is the amount of processed tasks of depth i
    atomic_llong tasksProcessed[maxDepth];
    // tasksToProcess[i] is the amount of tasks added but not yet processed
    atomic_llong tasksToProcess[maxDepth];

    // time point of the start of the program
    chrono::steady_clock::time_point start = chrono::steady_clock::now();

    // time point of last save
    chrono::steady_clock::time_point lastSave = chrono::steady_clock::now();

    // lastTaskTaken[i] is the time point of the latest moment when
    // a task of depth i has been taken
    chrono::steady_clock::time_point lastTaskTaken[maxDepth];
    mutex timeMutex;

    // return the amount of seconds elapsed since the timePoint
    int getSecondsElapsed(const chrono::steady_clock::time_point & timePoint) {
        return (int) chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - timePoint).count();
    }

    // print to the terminal the time elapsed since the timePoint
    void eprintTimeElapsed(const chrono::steady_clock::time_point & timePoint) {
        int seconds = getSecondsElapsed(timePoint);
        int minutes = seconds / 60, hours = minutes / 60;
        seconds %= 60, minutes %= 60;
        eprintf("%4d:%02d:%02d", hours, minutes, seconds);
    }

    // save to the file the amount of tasks processed and current stack for each depth
    void save() {
        lastSave = chrono::steady_clock::now();
        ofstream saveFile;
        saveFile.open("save.txt", ios_base::out);
        saveFile << SetPrecision(18, 15);
        for (int i = 0; i < maxDepth; ++i) {
            saveFile << tasksProcessed[i].load() << ' ';
        }
        saveFile << '\n';
        {
            unique_lock<mutex> lock(pool.tasksMutex);
            for (int i = sz(pool.tasks) - 1; i >= 0; --i) {
                for (auto t : pool.tasks[i]) {
                    saveFile << t << ' ' << i << '\n';
                }
            }
        }
        saveFile.close();
    }

    // once a second show current progress in the terminal
    void show() {
        while (true) {
            eprintf("\33[2J\33[H");
            eprintf("time elapsed: ");
            eprintTimeElapsed(start);
            eprintf("\nlast save was made %d second(s) ago\n", getSecondsElapsed(lastSave));

            bool noTasksLeft = true;

            {
                unique_lock<mutex> lock(timeMutex);
                for (int i = 0; i < maxDepth; ++i) {
                    if (tasksProcessed[i].load() != 0 || tasksToProcess[i].load() != 0) {
                        eprintf("depth %2d: %18lld task(s) processed, %18lld task(s) pending", i, tasksProcessed[i].load(), tasksToProcess[i].load());
                        if (tasksToProcess[i].load() != 0) {
                            noTasksLeft = false;
                            eprintf(", time since last task was taken:");
                            eprintTimeElapsed(lastTaskTaken[i]);
                            eprintf("\n");
                        } else {
                            eprintf("\n");
                        }
                    }
                }
            }

            // once a minute save current progress to the file
            if (getSecondsElapsed(lastSave) >= 60) {
                save();
            }

            if (noTasksLeft) {
                break;
            }
            this_thread::sleep_for(chrono::milliseconds(1000));
        }

        eprintf("\n");
    }

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
        cyclePart(const interval & x, const interval & y, const interval & alpha, Verifier *verifier) {
            interval beta = verifier->sumAlphaBeta - alpha;

            // W_1 and W_2
            Z_1W_1.o = point(x + cos(alpha), sin(alpha));
            Z_2W_2.o = point(sin(beta), y + cos(beta));
            
            Z_1W_1.v = point(cos(2 * alpha), sin(2 * alpha));
            Z_2W_2.v = point(sin(2 * beta), cos(2 * beta));
            V = intersect(Z_1W_1, Z_2W_2);

            len = verifier->d - Z_1W_1.o.x + verifier->d - Z_2W_2.o.y;
            len += (V - Z_1W_1.o).len() + (V - Z_2W_2.o).len();
        }
    };

    // initialize verifier
    void init(char** argv) {
        sscanf(argv[1], "%d", &threads);
        
        if (threads == 0) {
            // try to guess the number of concurrent threads supported by the machine
            threads = thread::hardware_concurrency() - 1;
        }
        
        // calculate L_0
        cyclePart optC(x_0, x_0, alpha_0, this);
        L_0 = optC.len + optC.V.len() - 1;

        // initialize delta and shift arrays
        deltas = vector<point6d>(maxDepth);
        deltas[0] = point6d(d - 1, d - 1, pi / 12, pi / 2, pi / 2, pi / 2) * cos60;
        for (int i = 1; i < maxDepth; ++i) {
            deltas[i] = deltas[i - 1] * cos60;
        }

        shift = vector<vector<point6d> >(maxDepth, vector<point6d>(64));
        for (int i = 0; i < maxDepth; ++i) {
            for (int j = 0; j < 64; ++j) {
                shift[i][j].x     = deltas[i + 1].x     * ((j & 1)  ? 1 : -1);
                shift[i][j].y     = deltas[i + 1].y     * ((j & 2)  ? 1 : -1);
                shift[i][j].alpha = deltas[i + 1].alpha * ((j & 4)  ? 1 : -1);
                shift[i][j].xi    = deltas[i + 1].xi    * ((j & 8)  ? 1 : -1);
                shift[i][j].xi_1  = deltas[i + 1].xi_1  * ((j & 16) ? 1 : -1);
                shift[i][j].xi_2  = deltas[i + 1].xi_2  * ((j & 32) ? 1 : -1);
            }
        }

        for (int i = 0; i < maxDepth; ++i) {
            lastTaskTaken[i] = start;
        }

        // initialize thread pool
        pool.init(threads, this);
    }

    // consider the box P(c, delta), where delta = deltas[depth];
    // if L(c) - Err(c, delta) < L_0, divide the box into 64 smaller boxes
    // and add a task for each of them;
    void checkBox(const point6d & c, Melzak & melzak, int depth) {
        point6d & delta = deltas[depth];

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
      
        cyclePart cycle(c.x, c.y, c.alpha, this);

        interval beta = sumAlphaBeta - c.alpha;

        point p[4];
        p[0] = cycle.V;
        p[1] = point(cos(c.xi), sin(c.xi));
        p[2] = point(c.x - cos(c.xi_1), sin(c.xi_1));
        p[3] = point(sin(c.xi_2), c.y - cos(c.xi_2));

        // length = L(c)
        interval length = melzak.solve(p) + cycle.len;

        // calculate Err(c, delta)
        interval M1 =  (4 * (c.x - delta.x) * cos(2 * (beta - delta.alpha)) + 4 * (c.y + delta.y) * sin(2 * (beta - delta.alpha)));
        M1 +=  (4 * sin(beta + delta.alpha) + 3 * cos(c.alpha + 2 * beta + delta.alpha) - cos(3 * c.alpha + 2 * beta + delta.alpha));
        interval M2 = -(4 * (c.x + delta.x) * cos(2 * (beta + delta.alpha)) + 4 * (c.y - delta.y) * sin(2 * (beta + delta.alpha)));
        M2 += -(4 * sin(beta - delta.alpha) + 3 * cos(c.alpha + 2 * beta - delta.alpha) - cos(3 * c.alpha + 2 * beta - delta.alpha));
        interval dV = mymax(interval(0), mymax(M1, M2)) * delta.alpha * interval(2) / interval(3);
        dV += delta.x * sin(2 * (c.alpha - delta.alpha)) / sin60 + delta.y * sin(2 * (beta - delta.alpha)) / sin60;
        interval Err = delta.xi_1 + delta.xi_2 + delta.xi + dV;

        Err += (1 - sin(2 * (c.alpha - delta.alpha) - 2 * pi / 3) / sin60) * delta.x;
        Err += (1 - cos(2 * (c.alpha + delta.alpha) - 2 * pi / 3) / sin60) * delta.y;

        interval F1 =  sin(c.alpha + delta.alpha) - sin(beta - delta.alpha);
        F1 += (2 * (c.x + delta.x) * cos(2 * (c.alpha - delta.alpha) - 2 * pi / 3) - 2 * (c.y - delta.y) * sin(2 * (c.alpha - delta.alpha) - 2 * pi / 3)) / sin60;
        F1 += (-sin(c.alpha - delta.alpha - pi / 6) + cos(c.alpha - delta.alpha - pi / 4)) / sin60;
        interval F2 = -sin(c.alpha - delta.alpha) + sin(beta + delta.alpha);
        F2 -= (2 * (c.x - delta.x) * cos(2 * (c.alpha + delta.alpha) - 2 * pi / 3) - 2 * (c.y + delta.y) * sin(2 * (c.alpha + delta.alpha) - 2 * pi / 3)) / sin60;
        F2 -= (-sin(c.alpha + delta.alpha - pi / 6) + cos(c.alpha + delta.alpha - pi / 4)) / sin60;
        Err += mymax(interval(0), mymax(F1, F2)) * delta.alpha;
        
        if (getSign(length - Err - L_0) == 1) {
            // the inequality is verified for every point in the box
            return;
        }

        // add new tasks
        tasksToProcess[depth + 1] += 64;
        for (int i = 0; i < 64; ++i) {
            point6d nc = c + shift[depth][i];
            pool.addTask(nc, depth + 1);
        }

        return;
    }

    // if there is no save file detected, add the initial task;
    // if there is a save file, load the amount of tasks processed and current stack for each depth
    void addTasks() {
        ifstream saveFile("save.txt");
        if (saveFile.is_open()) {
            string line;
            {
                getline(saveFile, line);
                istringstream iss(line);
                for (int i = 0; i < maxDepth; ++i) {
                    long long tmp;
                    iss >> tmp;
                    tasksProcessed[i] = tmp;
                }
            }
            while (getline(saveFile, line)) {
                istringstream iss(line);
                point6d t;
                int depth;
                iss >> t >> depth;
                pool.addTask(t, depth);
                ++tasksToProcess[depth];
            }
            saveFile.close();
        } else {
            point6d t = deltas[0];
            t.alpha += sumAlphaBeta - pi / 2;
            ++tasksToProcess[0];
            pool.addTask(t, 0);
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
    verifier.save();
    return 0;
}
