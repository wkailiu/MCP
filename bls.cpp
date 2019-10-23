#include <iostream>
#include <vector>
#include <algorithm>

#include <iterator>
#include <set>
#include <list>
#include <string>
#include <cstring>
#include <fstream>

#include <ctime>
#include <iomanip>

int N;
std::vector<float> W;

typedef std::vector<int> Clique;

class NodeSet {
	std::vector<bool> state;
	void init(int n) {
		state = std::vector<bool>(n, false);
	}
public:
	NodeSet(int n) { init(n); };
	NodeSet() { init(N); };
	NodeSet(Clique C) {
		init(N);
		for (auto v : C) {
			add(v);
		}
	}

	void add(int n) { state[n] = true; }
	void remove(int n) { state[n] = false; }
	bool have(int n) { return state[n]; }

	NodeSet set_intersection(NodeSet ns2) {
		NodeSet ns1;
		for (int i = 0; i < N; i++) {
			if (have(i) && ns2.have(i)) {
				ns1.add(i);
			}
		}
		return ns1;
	}

	void reverse() {
		for (int i = 0; i < N; i++) {
			state[i] = !state[i];
		}
	}

	NodeSet set_union(NodeSet ns2) {
		NodeSet ns1;
		for (int i = 0; i < N; i++) {
			if (have(i) || ns2.have(i)) {
				ns1.add(i);
			}
		}
		return ns1;
	}
	int size() {
		int cnt = 0;
		for (int i = 0; i < N; i++) {
			if (state[i])cnt++;
		}
		return cnt;
	}
	
};


class Edge {
public:
	int v;
	int u;
	Edge(int v0, int u0) {
		v = v0;
		u = u0; 
	}
};
typedef std::vector<Edge> EdgeSet;

class Graph {
public:
	int numEdge;
	std::vector<NodeSet> adjacentNode;
	std::vector<std::vector<bool>> isAdjacent;
	Clique MCP;
};

class Move {
public:
	int type;	//1:M1 2:M2 3:M3 4:M4
	int v;
	int u;
	Move(int type0, int v0, int u0 = -1) {
		type = type0;
		v = v0;
		u = u0;
	}
};
typedef std::vector<Move> MoveSet;

const float INF = 12345678.9f;
int Iter = 0;
Graph G;
std::vector<int> TL;		// tabu list
float f_best;
NodeSet PA(N);
EdgeSet OM;
NodeSet OC(N);
Clique C_best;
bool find_best_solution = false;


// 参数
int L_0;					// 0.01*N
const int T = 1000;
int L_max;				// 0.1*N
const float coefficient_strong = 0.8;
const int fa = 7;				//coefficient for tabu tenure
const float P0 = 0.75;
const float coefficient_random = 0.8;


void print(NodeSet ns) {
	std::cout << "[";
	for (int i = 0; i < N; i++) {
		if (ns.have(i)) {
			std::cout << i << " ";
		}
	}
	std::cout << "]" << std::endl;
}
void print(Graph g) {
	// 打印图的信息
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "Graph:" << std::endl;
	std::cout << "node:" << N << "\tedge:" << g.numEdge << std::endl;
	std::cout << "MCP[" << g.MCP.size() << "]:" << std::endl;
	for (int v : g.MCP) {
		std::cout << v << " ";
	}
	std::cout << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
}
void print(EdgeSet es) {
	std::cout << "[";
	for (auto e : es) {
		std::cout << "(" << e.v << "," << e.u << ") ";
	}
	std::cout << "]" << std::endl;
}
void print(Clique C) {
	std::cout << "Iter = " << Iter << "\tf_best=" << f_best << "\t[";
	for (auto c : C) {
		std::cout << c << " ";
	}
	std::cout << "]" << std::endl;
}

/*
input:		" 26 120 119 157 "
output:		[26, 120, 119, 157]
*/
std::vector<int> parseStringIntVector(std::string s) {
	std::vector<int> vec;
	int a = 0;
	bool start = false;
	for (char c:s) {
		if (c == ' ') {
			if (start) {
				vec.push_back(a);
				start = false;
				a = 0;
			}
		}
		else {
			a = a * 10 + c - '0';
			start = true;
		}
	}
	if (start) {
		vec.push_back(a);
	}
	return vec;
}

// read graph from *.cql file
Graph readClqFile(std::string filename) {
	Graph g;
	char head;
	std::string line;
	std::ifstream fr(filename, std::ios::in);
	while (!fr.eof()) {
		fr >> head;
		getline(fr, line);

		if (head == 'e') {
			std::vector<int> vec = parseStringIntVector(line);
			// 需要-1
			int u = vec[0] - 1;
			int v = vec[1] - 1;
			g.adjacentNode[u].add(v);
			g.adjacentNode[v].add(u);
			g.isAdjacent[u][v] = true;
			g.isAdjacent[v][u] = true;
		}

		else if (head == 'c') {
			if (line.find("Clique Elements are") != std::string::npos) {
				while (1) {
					fr >> head;
					getline(fr, line);
					if (line.find("Estimated Size") != std::string::npos)
						break;
					std::vector<int> vec = parseStringIntVector(line);
					for (int a : vec) {
						g.MCP.push_back(a);
					}
				}
			}
		}

		else if (head == 'p') {
			std::vector<int> vec = parseStringIntVector(line.substr(5));
			N = vec[0];
			g.numEdge = vec[1];
			for (int i = 0; i < N; i++) {
				g.adjacentNode.push_back(NodeSet(N));
				g.isAdjacent.push_back(std::vector<bool>(N, false));
			}
		}
	}

	return g;

}

// Calculate the weight of the graph
float f(Clique C) {
	return C.size();
	float w = 0;
	for (auto v : C) {
		w += W[v];
	}
	return w;
}

Clique generate_initial_solution() {
	Clique C;
	int v = rand() % N;
	C.push_back(v);
	while (1) {
		NodeSet candidate;
		for (int c : C) {
			candidate = candidate.set_union(G.adjacentNode[c]);
		}
		int n1 = C.size();
		for (int u = 0; u < N; u++) {
			if (!candidate.have(u))continue;
			bool isGood = true;
			for (int c : C) {
				if (!G.isAdjacent[u][c]) {
					isGood = false; 
					break;
				}
			}
			if (isGood) {
				C.push_back(u);
			}
		}

		if (C.size() == n1) {
			break;
		}
	}
	return C;
}

// update PA,OM,OC of the clique
void update_PA_OM_OC(Clique C) {
	//PA:自己不属于团,却与团内所有点相邻
	PA = NodeSet();
	for (int i = 0;i < N; i++) {
		bool isGood = true;
		for (auto v : C) {
			if (i == v || (!G.isAdjacent[i][v])) {
				isGood = false;
				break;
			}
		}
		if (isGood) {
			PA.add(i);
		}
	}
	/*
	OM: pairs(v,u)
		v not in C
		u in C
		v,u is not adjacent
		v与团内所有点相连除了u
	*/
	OM = EdgeSet();
	for (int v = 0; v < N; v++) {
		bool isGood = true;
		int goodU = -1;
		int cnt = 0;
		for (auto u : C) {
			if (v == u) {
				isGood = false; break;
			}
			if (!G.isAdjacent[v][u]) {
				cnt++;
				goodU = u;
				if (cnt > 1) {
					isGood = false; break;
				}
			}
		}
		if (isGood && goodU >= 0) {
			OM.push_back(Edge(v, goodU));
		}
	}

	// OC: 所有团外的点
	OC = NodeSet();
	for (auto v : C) {
		OC.add(v);
	}
	OC.reverse();
}

Clique move(Clique C, Move m) {
	Clique C2;
	switch (m.type)
	{
	case 1: {
		// add v
		for (auto v : C) {
			C2.push_back(v);
		}
		C2.push_back(m.v);
	}break;
	case 2: {
		// add v, remove u
		for (auto v : C) {
			if (v != m.u)
				C2.push_back(v);
		}
		C2.push_back(m.v);
	}break;
	case 3: {
		// remove v
		for (auto v : C) {
			if (v != m.v)
				C2.push_back(v);
		}
	}break;
	case 4: {
		// add v, repair C to clique
		for (auto v : C) {
			if (G.isAdjacent[v][m.v])
				C2.push_back(v);
		}
		C2.push_back(m.v);
	}break;
	default:
		break;
	}
	return C2;
}

bool prohibited(Clique C, Move m) {
	if (m.type == 1)return false;

	int gamma = fa;
	if (OM.size() > 0)
		gamma += rand() % OM.size();
	if (m.type == 4) {
		for (auto v : C) {
			if (!G.isAdjacent[v][m.v]) {
				if (TL[v] < 0) return false;
				return (Iter - TL[v]) <= gamma;
			}
		}
	}
	else {
		int v;
		if (m.type == 2) v = m.u;
		if (m.type == 3) v = m.v;
		if (TL[v] < 0) return false;
		return (Iter - TL[v]) <= gamma;
	}
}

MoveSet M1_M2(Clique C) {
	MoveSet M12;
	// M1: add v
	for (int v = 0; v < N; v++) {
		if (PA.have(v)) {
			M12.push_back(Move(1, v));
		}
	}
	// M2:
	for (auto e : OM) {
		M12.push_back(Move(2, e.v, e.u));
	}
	
	// Select the best move m from the set of moves formed by the union M1 ∪ M2
	MoveSet bestMove;
	float best_fc = -INF;
	for (Move m : M12) {
		Clique C2 = move(C, m);
		float fc = f(C2);
		if (fc > best_fc) {
			best_fc = fc;
			bestMove.clear();
			bestMove.push_back(m);
		}
		else if (fc == best_fc) {
			bestMove.push_back(m);
		}
	}
	if (best_fc <= f(C))
		bestMove.clear();
	return bestMove;
}

MoveSet M4_Set(Clique C,float alpha) {
	MoveSet M4;
	for (int v = 0; v < N; v++) {
		if (OC.have(v)) {
			float sum_w = W[v];
			for (auto u : C) {
				if (G.isAdjacent[v][u]) {
					sum_w += W[u];
				}
			}
			if (sum_w >= alpha*f(C)) {
				M4.push_back(Move(4, v));
			}
		}
	}
	return M4;
}

void update_tabu_list(Clique C_old, Clique C_new) {
	NodeSet ns_old(C_old);
	NodeSet ns_new(C_new);

	for (int i = 0; i < N; i++) {
		if (ns_old.have(i) && (!ns_new.have(i))) {
			TL[i] = Iter;
		}
	}
}

MoveSet A_Set(Clique C) {
	MoveSet A;
	// M1: add v
	for (int v = 0; v < N; v++) {
		if (PA.have(v)) {
			A.push_back(Move(1, v));
		}
	}
	// M2:
	for (auto e : OM) {
		A.push_back(Move(2, e.v, e.u));
	}
	// M3:
	for (auto v : C) {
		A.push_back(Move(3, v));
	}

	MoveSet bestMove;
	float best_fc = -INF;
	for (auto m : A) {
		Clique C2 = move(C, m);
		float best_fc = -INF;
		float fc = f(C2);
		if ((!prohibited(C, m)) || fc > f_best) {
			if (fc > best_fc) {
				best_fc = fc;
				bestMove.clear();
				bestMove.push_back(m);
			}
			else if (fc == best_fc){
				bestMove.push_back(m);
			}
		}
	}
	return bestMove;
}

Clique Perturb(Clique C, int L, MoveSet M) {
	if (M.size() == 0)return C;
	for (int i = 0; i < L; i++) {
		int k = rand() % M.size();
		Clique C2 = move(C, M[k]);
		update_tabu_list(C, C2);
		C = C2;
		float fc = f(C);
		if (fc > f_best) {
			f_best = fc;
			
			C_best = C;
			sort(C_best.begin(), C_best.end());
			print(C_best);
			if (C_best == G.MCP)
				find_best_solution = true;
		}
		update_PA_OM_OC(C);
		Iter += 1;
	}
	return C;
}

float Formula_1(int w) {
	float P = exp(-w*1.0 / T);
	return std::max(P, P0);
}

Clique Perturbation(Clique C, int L, int w) {
	if (w == 0) {
		//Strong random perturb
		//std::cout << "Strong random perturb" << std::endl;
		MoveSet M4 = M4_Set(C, coefficient_strong);
		return Perturb(C, L, M4);
	}
	else {
		float P = Formula_1(w);
		float p = (rand() % 10000)*1.0 / 10000;
		if (p <= P) {
			//Directed perturbation with moves from set A
			MoveSet A = A_Set(C);
			return Perturb(C, L, A);
		}
		else {
			//random perturb
			MoveSet M4 = M4_Set(C, coefficient_random);
			return Perturb(C, L, M4);
		}
	}
}


int main() {
	//srand((unsigned)time(NULL));
	std::string filename = "brock200_2.clq";	//1436ms
	G = readClqFile(filename);	
	
	std::cout << filename << std::endl;
	
	sort(G.MCP.begin(), G.MCP.end());
	print(G);

	clock_t t1 = clock();

	W = std::vector<float>(N, 1.0f);
	TL = std::vector<int>(N, -1);
	L_0 = (int)(0.01*N);
	L_max = (int)(0.1*N);

	Clique C = generate_initial_solution();
	
	update_PA_OM_OC(C);

	int L = L_0;
	float fc = f(C);
	C_best = C;
	f_best = fc;

	std::cout << "initial C: " << std::endl;
	print(C);
	std::cout << "PA: "; print(PA);
	std::cout << "OM: "; print(OM);
	std::cout << "OC: "; print(OC);

	Clique C_p = C;
	int w = 0;

	while (!find_best_solution) {
		while (true) {
			MoveSet M12 = M1_M2(C);
			if (M12.size() == 0)break;
			int k = rand() % M12.size();
			Move m = M12[k];

			Clique C2 = move(C, m);
			if (f(C2) > f(C)) {
				update_tabu_list(C, C2);
				C = C2;
				update_PA_OM_OC(C);
				Iter += 1;
			}
			else {
				break;
			}
		}
		fc = f(C);
		if (fc > f_best) {
			C_best = C; f_best = fc;

			sort(C_best.begin(), C_best.end());
			print(C_best);
			if (C_best == G.MCP)find_best_solution = true;
			w = 0;
		}
		else {
			w += 1;
		}
		if (w > T) {
			L = L_max;
			w = 0;
		}
		else if (C == C_p) {
			L += 1;
		}
		else {
			L = L_0;
		}
		C_p = C;
		C = Perturbation(C, L, w);
	}

	clock_t t2 = clock();
	std::cout << "time : "<<(t2 - t1) * 1.0 / CLOCKS_PER_SEC * 1000<<"ms" << std::endl;
	return 0;
}