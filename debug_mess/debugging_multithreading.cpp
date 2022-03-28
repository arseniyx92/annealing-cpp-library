#define ASYNC 1
#if ASYNC
#include <future>
#include <mutex>
static std::mutex Tmutex;
#else
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <cmath>
#include <thread>
using namespace std;

template <typename T>
struct random_generator {
    int type = 0; double arg;
    mt19937 rnd;
    uniform_real_distribution<double> probability_distribution;
    uniform_real_distribution<double> double_distribution;
    uniform_int_distribution<int> integer_distribution;
    // TYPES:
    //      0) usual mt19937 (normal distribution)
    //      1) BOLTZMANN ANNEALING (normal distribution)
    //      2) FAST ANNEALING (Cauchy distribution)
    //      3) uniform_real_distribution
    //      4) uniform_int_distribution
    random_generator(): type(0), arg(1) {
        probability_distribution = uniform_real_distribution<double>(0., 1.);
        //rnd.seed(chrono::steady_clock::now().time_since_epoch().count());
        rnd.seed(55);
    }
    random_generator(int type): type(type), arg(1) {
        probability_distribution = uniform_real_distribution<double>(0., 1.);
        //rnd.seed(chrono::steady_clock::now().time_since_epoch().count());
        rnd.seed(55);
    }
    random_generator(int type, double arg): type(type), arg(arg) {
        probability_distribution = uniform_real_distribution<double>(0., 1.);
        //rnd.seed(chrono::steady_clock::now().time_since_epoch().count());
        rnd.seed(55);
    }
    T gen(T x, double t, T l, T r) {
        T val;
        if (type == 0) {
            val = rnd();
        } else if (type == 1) {
            normal_distribution<> boltzmann(x, t*arg);
            val = boltzmann(rnd);
        } else if (type == 2) {
            cauchy_distribution<> cauchy(x, t*arg);
            val = cauchy(rnd);
        } else if (type == 3) {
            double_distribution = uniform_real_distribution<double>(l, r);
            val = double_distribution(rnd);
        } else if (type == 4) {
            integer_distribution = uniform_int_distribution<int>(l, r);
            val = integer_distribution(rnd);
        }
        if (l > val || r < val) {
            if (l >= 0 && val < 0) val *= -1.;
            val = val - (val / (r - l + 1)) * (r - l + 1);
            val += l;
        }
        return val;
    }
    double generate_probability() {
        return probability_distribution(rnd);
    }
}; // generator

struct decrease_func {
    int type = 2, k = 0; double d = 0.99, t0 = 0.9;
    // descent type
    //      0) BOLTZMANN DESCENT T[k] = T[0]/ln(k+1)
    //      1) CAUCHY DESCENT T[k] = T[0]/k
    //      2) EXPONENTIAL DESCENT T[k] = T[k-1]*d; where d [0;1] is a parameter
    // d parameter (optional if EXPONENTIAL DESCENT chosen)
    decrease_func(): k(0), d(0.99) {}
    decrease_func(double t0): k(0), d(0.99), t0(t0) {}
    decrease_func(int type, double t0): type(type), k(0), d(0.99), t0(t0) {}
    decrease_func(int type, double d, double t0): type(type), k(0), d(d), t0(t0) {}
    double decrease(double t) {
        k++;
        if (type == 0) {
            return t0/log(k+1);
        } else if (type == 1) {
            return t0/k;
        } else if (type == 2) {
            return t*d;
        }
        cerr << "INAPPROPRIATE DECREASE FUNCTION TYPE VALUE, should be 0 or 1 or 2" << endl;
        return t;
    }
}; // descent

struct acceptance_func {
    int type = 1; double eps = 1;
    // acceptance type (h)
    //      0) h = exp(-eps*(E[t]-E[t-1])/T) where eps is a parameter
    //      1) h = T
    // eps parameter (optional if h = exp(-e*(E[t]-E[t-1])/T) chosen)
    acceptance_func(): type(1) {}
    acceptance_func(int type): type(type), eps(1) {}
    acceptance_func(int type, double eps): type(type), eps(eps) {}
    bool decide(double delta, double t, double x) {
        double h = 0.;
        if (type == 0) {
            h = exp(-delta*eps/t);
        } else if (type == 1) {
            h = t;
        }
        if (x < h) {
            return true;
        } else {
            return false;
        }
    }
}; // state acceptance

template <typename T, typename G>
struct State {
    G f;
    int n;
    vector<int> perm;
    State() = default;
    State(int n): n(n) {
        mt19937 rnd(55);
        perm.resize(n); iota(perm.begin(), perm.end(), 0); shuffle(perm.begin(), perm.end(), rnd);
        f = F();
    };
    State(int n, vector<int> perm): n(n), perm(move(perm)) {
        f = F();
    }
    ~State() {
        perm.clear();
        perm.shrink_to_fit();
    }
    State<T, G> generate_new_state(random_generator<T>& gen, double t) {
        vector<int> nperm = perm;
        int i = gen.gen(n/2, t, 0, n-1);
        int j = gen.gen(n/2, t, 0, n-1);
        swap(nperm[i], nperm[j]);
        return State<T, G>(n, nperm);
    }
    T F() {
        f = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = i+1; j < n; ++j) {
                if ((j-i)+perm[i] == perm[j] || (i-j)+perm[i] == perm[j]) f++;
            }
        }
        return f;
    }
}; // TODO this structure needs user's alterations

template <typename T, typename G>
struct annealizer { // only for minimization (if you need to maximize change all functions out from x to -x)
    int N; double t;
    State<T, G> best_state, current_state;
    random_generator<T> gen;
    decrease_func descent;
    acceptance_func acceptance;
    annealizer(State<T, G> initial_state, int iterations, double temperature, int generator_type, double generator_arg, int descent_type, double descent_arg, int AC_type, double AC_arg) {
        N = iterations;
        t = temperature;
        gen = random_generator<T>(generator_type, generator_arg);
        descent = decrease_func(descent_type, descent_arg, t);
        acceptance = acceptance_func(AC_type, AC_arg);
        best_state = current_state = initial_state;
    }
    State<T, G> anneal() {
        for (int iter = 0; iter < N; ++iter) {
            t = descent.decrease(t);
            State<T, G> new_state = current_state.generate_new_state(gen, t);
            if ((new_state.f < current_state.f) || (acceptance.decide(new_state.f-current_state.f, t, gen.generate_probability()))) {
                current_state = new_state;
                if (new_state.f < best_state.f) best_state = new_state;
            }
        }
        return best_state;
    };
    //1 initial guess
    //1 estimation function F()
    //1 recalculating function
    //1 alteration function
    //1 generator function
    //1 number of iterations
    //1 initial temperature
    //1 generator type
    //      0) pure randomness
    //      1) BOLTZMANN ANNEALING
    //      2) FAST ANNEALING (Cauchy)
    //1 descent type
    //      0) BOLTZMANN DESCENT T[k] = T[0]/ln(k+1)
    //      1) CAUCHY DESCENT T[k] = T[0]/k
    //      2) EXPONENTIAL DESCENT T[k] = T[k-1]*d; where d [0;1] is a parameter
    //1 d parameter (optional if EXPONENTIAL DESCENT chosen)
    //1 acceptance type (h)
    //      0) h = exp(-eps*(E[t]-E[t-1])/T) where eps is a parameter
    //      1) h = T
    //1 eps parameter (optional if h = exp(-e*(E[t]-E[t-1])/T) chosen)
    //0 savings file path
};

template <typename T, typename G>
static void fixedTemperatureForSearch(State<T, G> initial_state, int iterations, double temperature, double accept_arg, vector<pair<double, double> >* TemperatureResults) {
    annealizer instance(initial_state, iterations, temperature, 4, 1., 2, 0.9, 0, accept_arg);
#if ASYNC
    lock_guard<mutex> lock(Tmutex);
#else
#endif
    TemperatureResults->emplace_back(instance.anneal().f, temperature);
    //cerr << this_thread::get_id() << endl;
}

template <typename T, typename G>
static void fixedDescentForSearch(State<T, G> initial_state, int iterations, double temperature, double accept_arg, double descent, vector<pair<double, double> >* DescentResults) {
    annealizer instance(initial_state, iterations, temperature, 4, 1., 2, descent, 0, accept_arg);
#if ASYNC
    lock_guard<mutex> lock(Tmutex);
#else
#endif
    DescentResults->emplace_back(instance.anneal().f, descent);
}

static void keks() {
    int x = 0;
    x++;
}

template <typename T, typename G>
State<T, G> autoSearch(State<T, G> initial_state, double SecondsToWait = 30.) {
    // achieving appropriate number of iterations
    double TimeForEachIteration = clock();
    random_generator<T> SampleGenerator(4);
    initial_state.generate_new_state(SampleGenerator, 1.);
    TimeForEachIteration = ((double)clock()-TimeForEachIteration)/CLOCKS_PER_SEC;
    const double localTime = SecondsToWait/5.;
    int APPROBATION_ITERATIONS = ceil(localTime/TimeForEachIteration);
    int ITERATIONS = ceil(SecondsToWait/TimeForEachIteration);
    cerr << "FOUND BEST ITERATIONS: " << ITERATIONS << endl;

    // achieving appropriate accept function
    double ACCEPT_ARG = 2; int cur = (int)initial_state.f;
    while (cur > 0) {
        cur /= 10;
        ACCEPT_ARG /= 10.;
    }
    cerr << "FIXED ACCEPT TYPE: " << 0 << endl;
    cerr << "FOUND BEST ACCEPT ARGUMENT: " << ACCEPT_ARG << endl;

    // achieving appropriate temperature
    random_generator<double> TemperatureGenerator(3);
    const int K = 50;
#if ASYNC
    vector<future<void> > TemperatureFutures;
#else
#endif
    vector<pair<double, double> > TemperatureResults;
    for (int i = 0; i < K; ++i) {
        double t = TemperatureGenerator.gen(0., 1., 0.65, 0.95);
#if ASYNC
        TemperatureFutures.push_back(async(launch::async, fixedTemperatureForSearch<T, G>, initial_state, APPROBATION_ITERATIONS, t, ACCEPT_ARG, &TemperatureResults));
#else
        fixedTemperatureForSearch(initial_state, APPROBATION_ITERATIONS, t, ACCEPT_ARG, &TemperatureResults);
#endif
    }
#if ASYNC
//    std::this_thread::sleep_for(std::chrono::seconds(1000));
//    cout << "KEK" << endl;
    for (int i = 0; i < K; ++i) {
        TemperatureFutures[i].wait();
    }
#else
#endif
    sort(TemperatureResults.begin(), TemperatureResults.end());
    double TEMPERATURE = TemperatureResults[0].second;
    cerr << "FOUND BEST TEMPERATURE: " << TEMPERATURE << endl;

    // achieving appropriate descent function
    random_generator<double> DescentGenerator(3);
    const int H = 50;
#if ASYNC
    vector<future<void> > DescentFutures;
#else
#endif
    vector<pair<double, double> > DescentResults;
    for (int i = 0; i < H; ++i) {
        double descent = DescentGenerator.gen(0., 1., 0.8, 0.999);
#if ASYNC
        DescentFutures.push_back(async(launch::async, fixedDescentForSearch<T, G>, initial_state, APPROBATION_ITERATIONS, TEMPERATURE, ACCEPT_ARG, descent, &DescentResults));
#else
        fixedDescentForSearch(initial_state, APPROBATION_ITERATIONS, TEMPERATURE, ACCEPT_ARG, descent, &DescentResults);
#endif
    }
#if ASYNC
    for (int i = 0; i < H; ++i) DescentFutures[i].wait();
#else
#endif
    sort(DescentResults.begin(), DescentResults.end());
    double DESCENT_ARG = DescentResults[0].second;
    cerr << "FIXED DESCENT TYPE: " << 2 << endl;
    cerr << "FOUND BEST DESCENT ARGUMENT: " << DESCENT_ARG << endl;

    // achieving appropriate generator arg
    // achieving appropriate generator type
    cerr << "FIXED GENERATOR TYPE: " << 4 << endl;
    cerr << "FIXED GENERATOR ARGUMENT: " << 1 << endl;

    cerr << "(State<T, G>, " << ITERATIONS << ", " << TEMPERATURE << ", 4, 1., 2, " << DESCENT_ARG << ", 0, " << ACCEPT_ARG << ")" << endl;
    annealizer instance(initial_state, ITERATIONS, TEMPERATURE, 4, 1., 2, DESCENT_ARG, 0, ACCEPT_ARG);
    return instance.anneal();
}


int main() {
    int n;
    cin >> n;
    State<int, int> instance(n);
    instance = autoSearch(instance, 2.);
    //annealizer<int, int> descent(instance, 30000, 0.9, 4, 1, 2, 0.99, 0, 0.02);
    //instance = descent.anneal();
    cout << instance.f << '\n';
    for (int x : instance.perm) cout << x + 1 << ' ';
    cout.precision(9);
    cout << endl << (double) clock() / CLOCKS_PER_SEC << endl;
}
// 100
// NO ASYNC 11.596006 sec
// WITH ASYNC 10.136434 sec

// 200
// NO ASYNC 16.726927 sec
// WITH ASYNC 14.547082 sec

// 500
// NO ASYNC 36.203969 sec
// ASYNC 30.946662 sec