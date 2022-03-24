#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <cmath>
using namespace std;

template <typename T>
struct random_generator {
    int type = 0;
    mt19937 rnd;
    uniform_real_distribution<double> probability_distribution;
    // TYPES:
    //      0) usual mt19937 (normal distribution)
    //      1) BOLTZMANN ANNEALING (normal distribution)
    //      2) FAST ANNEALING (Cauchy distribution)
    random_generator() {
        probability_distribution = uniform_real_distribution<double>(0., 1.);
        rnd.seed(chrono::steady_clock::now().time_since_epoch().count());
    }
    random_generator(int type): type(type) {
        probability_distribution = uniform_real_distribution<double>(0., 1.);
        rnd.seed(chrono::steady_clock::now().time_since_epoch().count());
    }
    T gen(T x, T t, T l, T r) {
        T val;
        if (type == 0) {
            val = rnd();
        } else if (type == 1) {
            normal_distribution<T> boltzmann(x, t);
            val = boltzmann(rnd);
        } else if (type == 2) {
            cauchy_distribution<T> cauchy(x, t);
            val = cauchy(rnd);
        }
        if (l > val) val = l;
        if (r < val) val = r;
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
    random_generator<T> gen;
    // TODO data maintainer
    //..
    //
    State() {
        // TODO State constructor
        //..
        //
    };
    State<T, G> generate_new_state() {
        // TODO new state generation
        //..
        //
    }
    T F() {
        // TODO new estimation function
        //..
        //
    }
}; // TODO this structure needs user's alterations

template <typename T, typename G>
struct annealizer { // only for minimization (if you need to maximize change all functions out from x to -x)
    int N; double t;
    State<T, G> best_state, current_state;
    random_generator<T> gen;
    decrease_func descent;
    acceptance_func acceptance;
    annealizer(State<T, G> initial_state, int iterations, double temperature, int generator_type, int descent_type, double descent_arg, int AC_type, double AC_arg) {
        N = iterations;
        t = temperature;
        gen = random_generator<T>(generator_type);
        descent = decrease_func(descent_type, descent_arg, t);
        acceptance = acceptance_func(AC_type, AC_arg);
        G init = F(initial_state);
        best_state = current_state = {init, initial_state};
    }
    State<T, G> anneal() {
        for (int iter = 0; iter < N; ++iter) {
            t = descent.decrease(t);
            State<T, G> new_state = current_state.generate_new_state(gen);
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
    //      1) BOLTZMANN ANNEALING
    //      2) FAST ANNEALING (Cauchy)
    //      3) pure randomness
    //1 descent type
    //      1) BOLTZMANN DESCENT T[k] = T[0]/ln(k+1)
    //      2) CAUCHY DESCENT T[k] = T[0]/k
    //      3) EXPONENTIAL DESCENT T[k] = T[k-1]*d; where d [0;1] is a parameter
    //1 d parameter (optional if EXPONENTIAL DESCENT chosen)
    //1 acceptance type (h)
    //      1) h = exp(-eps*(E[t]-E[t-1])/T) where eps is a parameter
    //      2) h = distrib(rnd(), 0., 1.)
    //1 eps parameter (optional if h = exp(-e*(E[t]-E[t-1])/T) chosen)
    //0 savings file path
};