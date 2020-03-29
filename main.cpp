// #define NDEBUG // uncomment to turn off assert checks (increases speed)

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
#include <iostream>
#include <string>
#include <set>

#define eprintf(...) fprintf(stderr, __VA_ARGS__), fflush(stderr)
#define sz(a) ((int) (a).size())

const int maxDepth = 20;
const int maxThreads = 64;

using namespace std;

// usage: ./main A_1Degs k cutoffDepth radiusXY radiusAlpha

// the program verifies that for every point (x, y, alpha, xi_up, xi_1, xi_2) in the parametric space,
// such that |x - x_0| > radiusXY, |y - y_0| > radiusXY, |alpha - alpha_0| > radiusAlpha,
// H(Sigma \cap ANGLE) >= H(Sigma_0 \cap ANGLE)
//
// A_1Degs is angle A_1 in degrees
//
// the parametric space is initially divided into k^6 hyperrectangles
// each of these initial hyperrectangles is considered as a separate task
// tasks are run in parallel using multiple threads
//
// cutoffDepth is the maximum allowed recursive depth
// (if this depth was reached by any hyperrectangle, the program needs to be
// restarted with higher cutoffDepth)

// progress bar format:
// time|fin/all|thread|thread|...|thread|0:cnt[0] 1:cnt[1] ... d:cnt[d]
//
// time is the time elapsed since the start of the program
//
// fin is the number of finished tasks
// all is the total number of tasks
//
// |thread| is the information about one of the threads in format |d time|,
// where d is the current recursive depth in this thread,
// time is the time elapsed since the start of the last task in this thread
//
// cnt[i] is the number of finished tasks which had maximum reached recursive depth i

// a point in the parametric space
struct point6d {
    long double x, y, alpha, xi_up, xi_1, xi_2;

    point6d() {}
    point6d(const long double & x, const long double & y, const long double & alpha,
            const long double & xi_up, const long double & xi_1, const long double & xi_2):
                x(x), y(y), alpha(alpha), xi_up(xi_up), xi_1(xi_1), xi_2(xi_2) {}

    // return the point, multiplied by the constant
    point6d operator *(const long double & m) const {
        return point6d(x * m, y * m, alpha * m, xi_up * m, xi_1 * m, xi_2 * m);
    }

    // return the sum of two points
    point6d operator +(const point6d & p) const {
        return point6d(x + p.x, y + p.y, alpha + p.alpha, xi_up + p.xi_up, xi_1 + p.xi_1, xi_2 + p.xi_2);
    }
};

struct Verifier {
    // angle A_1 in degrees and radians
    long double A_1Degs, A_1;

    // direction vectors of lines A_1A_2 and A_1A_3
    point A_1A_2Dir, A_1A_3Dir;
    
    long double radiusXY, radiusAlpha;

    // sum of alpha and beta
    long double sumAlphaBeta;

    // d = |A_1D_11| = |A_1D_12|
    long double d;

    int k, cutoffDepth;

    // H(Sigma_0 \cap ANGLE), x_0, y_0, and alpha_0
    long double minLength, x_0, y_0, alpha_0;

    // deltas[i] stores the delta values for the hyperrectangle,
    // which is 2^i times smaller in each dimension than the initial hyperrectangle
    vector<point6d> deltas;

    // shift[i] is the list of the vectors from the center of the hyperrectangle,
    // which is 2^i times smaller in each dimension than the initial hyperrectangle,
    // to the centers of 64 smaller hyperrectangles into which it is divided
    vector<vector<point6d> > shift;

    // thread pool maintains the queue of tasks and distributes them among threads
    struct ThreadPool {
        vector<thread> threadPool;
        queue<point6d> tasks;
        mutex tasksMutex;
        condition_variable condition;
        atomic_bool terminate;
        atomic_bool stopped;
        
        void init(int threads, Verifier *verifier) {
            terminate = false;
            stopped = false;

            for (int i = 0; i < threads; i++) {
                threadPool.emplace_back(thread(&ThreadPool::Runner, this, i, verifier));
            }
        }

        ~ThreadPool() {
            if (!stopped) {
                ShutDown();
            }
        }

        // add new task to queue
        void Enqueue(point6d task) {
            {
                unique_lock<mutex> lock(tasksMutex);
                tasks.push(task);
            }

            // wake up one (random) thread for task.
            condition.notify_one();
        }

        void ShutDown() {
            terminate = true;
            condition.notify_all();

            for (thread &thread : threadPool) {
                thread.join();
            }

            threadPool.clear();
            stopped = true;
        }

        void Runner(int id, Verifier *verifier) {
            verifier->progressBar.currentDepth[id] = -1;
            point6d task;

            while (true) {
                {
                    unique_lock<mutex> lock(tasksMutex);

                    condition.wait(lock, [this] { return !tasks.empty() || terminate; });

                    if (terminate && tasks.empty()) {
                        return;
                    }

                    task = tasks.front();
                    tasks.pop();
                }

                try {
                    verifier->verify(task, id);
                } catch (runtime_error &ex) {
                    cerr << "at " << this_thread::get_id() << " " << ex.what() << endl;
                }
            }
        }
    } pool;

    // shows progress bar in the terminal
    struct ProgressBar {
        // the number of threads
        int threads;

        // currentDepth[i] is the current depth of recursion in thread i
        atomic_int currentDepth[maxThreads];

        // the amount of finished tasks
        atomic_int tasksFinishedAnyDepth;

        // tasksFinished[i] is the number of finished tasks which had maximum reached recursive depth i
        atomic_int tasksFinished[maxDepth];

        // total amount of tasks
        int tasksTotal;

        chrono::steady_clock::time_point start = chrono::steady_clock::now();
        mutex timeMutex;
        chrono::steady_clock::time_point threadStart[maxThreads];

        // show the time elapsed since the start of the program (if id = -1),
        // or since the start of the last task in thread id (if id >= 0)
        void eprintTimeElapsed(int id) {
            unique_lock<mutex> lock(timeMutex);

            int seconds = (int) chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - (id == -1 ? start : threadStart[id])).count();    
            int minutes = seconds / 60, hours = minutes / 60;
            seconds %= 60, minutes %= 60;
            if (hours > 0) {
                eprintf("%02d:", hours);
            }
            eprintf("%02d:%02d", minutes, seconds);
        }

        // thread id started the new task
        void setThreadStart(int id) {
            unique_lock<mutex> lock(timeMutex);
            threadStart[id] = chrono::steady_clock::now();
        }

        // every 500 milliseconds show current progress in the terminal
        void show(Verifier *verifier) {
            while (true) {
                eprintf("\33[2K\r");
                eprintTimeElapsed(-1);
                eprintf("|%d/%d|", tasksFinishedAnyDepth.load(), tasksTotal);
                for (int i = 0; i < threads; ++i) {
                    if (currentDepth[i].load() == -1) {
                        eprintf("----------|");
                        continue;
                    }
                    eprintf("%d ", currentDepth[i].load());
                    eprintTimeElapsed(i);
                    eprintf("|");
                }
                for (int i = 0; i < 20; ++i) {
                    if (i == verifier->cutoffDepth) {
                        eprintf("\033[1;31m");
                    }
                    if (tasksFinished[i].load() != 0) {
                        eprintf("%d:%d ", i, tasksFinished[i].load());
                    }
                }
                eprintf("\033[0m");
                if (tasksFinishedAnyDepth.load() == tasksTotal) {
                    break;
                }
                this_thread::sleep_for(chrono::milliseconds(500));
            }

            eprintf("\n");
        }
    } progressBar;


    // the part of the cycle C: Z_1-W_1-V-W_2-Z_2
    struct cyclePart {
        // lines through segments Z_1W_1 and Z_2W_2;
        // Z_1W_1.o = W_1, Z_2W_2.o = W_2
        line Z_1W_1, Z_2W_2;
        point V;

        // total length of Z_1-W_1-V-W_2-Z_2
        long double len;

        cyclePart() {}

        // find lines Z_1W_1, Z_2W_2 and point V, given x, y, alpha
        cyclePart(const long double & x, const long double & y, const long double & alpha, Verifier *verifier) {
            long double beta = verifier->sumAlphaBeta - alpha;

            // W_1 and W_2
            Z_1W_1.o = verifier->A_1A_2Dir * x + verifier->A_1A_2Dir.rotate(alpha);
            Z_2W_2.o = verifier->A_1A_3Dir * y + verifier->A_1A_3Dir.rotate(-beta);
            
            Z_1W_1.v = verifier->A_1A_2Dir.rotate(2 * alpha);
            Z_2W_2.v = verifier->A_1A_3Dir.rotate(-2 * beta);
            V = intersect(Z_1W_1, Z_2W_2);

            len = verifier->d - (Z_1W_1.o * verifier->A_1A_2Dir);
            len += verifier->d - (Z_2W_2.o * verifier->A_1A_3Dir);
            len += (V - Z_1W_1.o).len() + (V - Z_2W_2.o).len();
        }
    };

    void init(char** argv) {
        sscanf(argv[1], "%Lf", &A_1Degs);
        sscanf(argv[2], "%d", &k);
        sscanf(argv[3], "%d", &cutoffDepth);
        sscanf(argv[4], "%Lf", &radiusXY);
        sscanf(argv[5], "%Lf", &radiusAlpha);
        
        A_1 = A_1Degs / 180 * pi;
        A_1A_2Dir = point(1, 0);
        A_1A_3Dir = point(cos(A_1), sin(A_1));
        sumAlphaBeta = pi * 2 / 3 + A_1 / 2;
        d = 1.0 / sin60 / sin(A_1 / 2) + 2 + cos(A_1 / 2) / sin(A_1 / 2) + 1;

        x_0 = y_0 = 1 / sin(A_1 / 2);
        alpha_0 = sumAlphaBeta / 2;
        cyclePart optC(x_0, x_0, alpha_0, this);
        minLength = optC.len + optC.V.len() - 1;

        // print H(Sigma_0 \cap ANGLE) to the terminal
        eprintf("H(Sigma_0 \\cap ANGLE): %Lf\n", minLength);

        deltas = vector<point6d>(cutoffDepth);
        deltas[0] = point6d(d - 1, d - 1, pi - sumAlphaBeta, A_1, pi / 2, pi / 2) * (1.0 / k / 2);
        for (int i = 1; i < cutoffDepth; ++i) {
            deltas[i] = deltas[i - 1] * 0.5;
        }

        shift = vector<vector<point6d> >(cutoffDepth - 1, vector<point6d>(64));
        for (int i = 0; i < cutoffDepth - 1; ++i) {
            for (int j = 0; j < 64; ++j) {
                shift[i][j].x     = deltas[i + 1].x     * ((j & 1)  ? 1 : -1);
                shift[i][j].y     = deltas[i + 1].y     * ((j & 2)  ? 1 : -1);
                shift[i][j].alpha = deltas[i + 1].alpha * ((j & 4)  ? 1 : -1);
                shift[i][j].xi_up = deltas[i + 1].xi_up * ((j & 8)  ? 1 : -1);
                shift[i][j].xi_1  = deltas[i + 1].xi_1  * ((j & 16) ? 1 : -1);
                shift[i][j].xi_2  = deltas[i + 1].xi_2  * ((j & 32) ? 1 : -1);
            }
        }

        progressBar.threads = thread::hardware_concurrency() - 1;
        pool.init(progressBar.threads, this);
    }

    // consider the hyperrectangle [t.x - delta.x, t.x + delta.x] x [t.y - delta.y, t.y + delta.y] x
    // [t.alpha - delta.alpha, t.alpha + delta.alpha] x [t.xi_up - delta.xi_up, t.xi_up + delta.xi_up] x
    // [t.xi_1 - delta.xi_1, t.xi_1 + delta.xi_1] x [t.xi_2 - delta.xi_2, t.xi_2 + delta.xi_2],
    // where delta = deltas[depth];
    // if L(t) - error >= minLength, divide the hyperrectangle into 64 smaller hyperrectangles and
    // call this function recursively for each of them
    //
    // return maximum reached depth
    //
    // stop if recursion depth reaches cutoffDepth
    int checkHyperrectangle(const point6d & t, Melzak & melzak, int depth, int id) {
        progressBar.currentDepth[id] = depth;

        point6d & delta = deltas[depth];

        if (abs(t.x - x_0) + delta.x < radiusXY && abs(t.y - y_0) + delta.y < radiusXY && abs(t.alpha - alpha_0) + delta.alpha < radiusAlpha) {
            // |x - x_0| < radiusXY, |y - y_0| < radiusXY, and |alpha - alpha_0| < radiusAlpha
            // for each point in the hyperrectangle
            return depth;
        }
        if ((t.x + delta.x) * (t.x + delta.x) + (t.y + delta.y) * (t.y + delta.y) - 2 * (t.x + delta.x) * (t.y + delta.y) * cos(A_1) < 3.99) {
            // B_r(y(W_1)) intersects B_r(y(W_2)) for each point in the hyperrectangle
            return depth;
        }
      
        // calculate L(t)
        cyclePart b(t.x, t.y, t.alpha, this);

        long double beta = sumAlphaBeta - t.alpha;

        point p[4];
        p[0] = b.V;
        p[1] = point(cos(t.xi_up), sin(t.xi_up));
        p[2] = A_1A_2Dir * t.x + point(cos(t.xi_1), sin(t.xi_1));
        p[3] = A_1A_3Dir * t.y + point(cos(t.xi_2), sin(t.xi_2));

        long double length = melzak.solve(p) + b.len;

        // calculate error
        long double error = delta.xi_1 + delta.xi_2 + delta.xi_up;
        long double dV = 2 / sin60 * ((b.Z_1W_1.o - b.Z_2W_2.o).len() + (cos(t.alpha - delta.alpha) + cos(beta - delta.alpha)) / 2) * delta.alpha;
        dV += delta.x * sin(2 * (t.alpha - delta.alpha)) / sin60 + delta.y * sin(2 * (beta - delta.alpha)) / sin60;
        error += dV;
        long double dW_1 = delta.x + delta.alpha, dW_2 = delta.y + delta.alpha;
        error += 2 * (dW_1 * cos(t.alpha - delta.alpha) + dW_2 * cos(beta - delta.alpha) + dV * cos60);

        if (length - error > minLength) {
            // the inequality is verified for every point in the hyperrectangle
            return depth;
        }

        if (depth == cutoffDepth - 1) {
            // cutoffDepth reached; stop
            return cutoffDepth;
        }

        int reachedDepth = depth;

        // make recursive calls
        for (int i = 0; i < 64; ++i) {
            point6d nt = t + shift[depth][i];

            reachedDepth = max(reachedDepth, checkHyperrectangle(nt, melzak, depth + 1, id));
            if (reachedDepth == cutoffDepth) {
                // cutoffDepth reached; stop
                return reachedDepth;
            }
        }

        return reachedDepth;
    }

    // verify the initial hyperrectangle with center at t;
    // id is the identifier of the used thread
    void verify(point6d t, int id) {
        Melzak melzak;
        progressBar.setThreadStart(id);
        int reachedDepth = checkHyperrectangle(t, melzak, 0, id);
        progressBar.currentDepth[id] = -1;
        ++progressBar.tasksFinished[reachedDepth];
        ++progressBar.tasksFinishedAnyDepth;
    }

    // create a task for each initial hyperrectangle;
    // thanks to symmetry it is enough to consider only hyperrectangles each point in which satisfies x >= y
    void addTasks() {
        point6d t;
        
        progressBar.tasksTotal = 0;
        for (int xi = 0; xi < k; ++xi) {
            t.x = deltas[0].x * (2 * xi + 1);
            for (int yi = 0; yi <= xi; ++yi) {
                t.y = deltas[0].y * (2 * yi + 1);
                for (int alphai = 0; alphai < k; ++alphai) {
                    t.alpha = sumAlphaBeta - pi / 2 + deltas[0].alpha * (2 * alphai + 1);
                    for (int xi_upi = 0; xi_upi < k; ++xi_upi) {
                        t.xi_up = deltas[0].xi_up * (2 * xi_upi + 1);
                        for (int xi_1i = 0; xi_1i < k; ++xi_1i) {
                            t.xi_1 = pi / 2 + deltas[0].xi_1 * (2 * xi_1i + 1);
                            for (int xi_2i = 0; xi_2i < k; ++xi_2i) {
                                t.xi_2 = pi + A_1 + deltas[0].xi_2 * (2 * xi_2i + 1);
                                ++progressBar.tasksTotal;
                                pool.Enqueue(t);
                            }
                        }
                    }
                }
            }
        }
    }
} verifier;


int main(int argc, char** argv) {
    if (argc != 6) {
        eprintf("Wrong number of arguments\n");
        exit(1);
    }

    verifier.init(argv);
    verifier.addTasks();
    verifier.progressBar.show(&verifier);
    return 0;
}
