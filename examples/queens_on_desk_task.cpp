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
    State<T, G> generate_new_state(random_generator<T>& gen, double t) {
        vector<int> nperm = perm;
        int i = gen.gen(n/2, t, 0, n);
        int j = gen.gen(n/2, t, 0, n);
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

int main() {
    int n; cin >> n;
    State<int, int> instance(n);
    annealizer<int, int> descent(instance, 30000, 0.9, 4, 1, 2, 0.99, 0, 0.02);
    instance = descent.anneal();
    cout << instance.f << '\n';
    for (int x : instance.perm) cout << x+1 << ' ';
}