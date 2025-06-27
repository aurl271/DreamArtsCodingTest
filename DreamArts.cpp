// 福島直之DreamArtsコーディングテスト
// 駅と駅を結ぶ鉄道路線網において、もっとも長い片道きっぷの経路
// 同じ点を2回通ることはできないが、特別に始点と終点を同じにできるので1->2->1のようなケースも問題ないと想定
// グラフを連結な成分ごとに分けて考える
// 連結グラフが木構造かを判定し、木構造ならO(N)のアルゴリズムで最長経路を求め、木構造でなければビットDPでO(2^N*N^2)で求める

#include <iostream>
#include <vector>
#include <tuple>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <cassert>
#include <stack>

using namespace std;

// 駅を連結成分で分けたいので、union_findを使用する。
// ac-libraryから拝借
struct dsu {
  public:
    dsu() : _n(0) {}
    dsu(int n) : _n(n), parent_or_size(n, -1) {}

    int merge(int a, int b) {
        assert(0 <= a && a < _n);
        assert(0 <= b && b < _n);
        int x = leader(a), y = leader(b);
        if (x == y) return x;
        if (-parent_or_size[x] < -parent_or_size[y]) std::swap(x, y);
        parent_or_size[x] += parent_or_size[y];
        parent_or_size[y] = x;
        return x;
    }

    bool same(int a, int b) {
        assert(0 <= a && a < _n);
        assert(0 <= b && b < _n);
        return leader(a) == leader(b);
    }

    int leader(int a) {
        assert(0 <= a && a < _n);
        if (parent_or_size[a] < 0) return a;
        return parent_or_size[a] = leader(parent_or_size[a]);
    }

    int size(int a) {
        assert(0 <= a && a < _n);
        return -parent_or_size[leader(a)];
    }

    std::vector<std::vector<int>> groups() {
        std::vector<int> leader_buf(_n), group_size(_n);
        for (int i = 0; i < _n; i++) {
            leader_buf[i] = leader(i);
            group_size[leader_buf[i]]++;
        }
        std::vector<std::vector<int>> result(_n);
        for (int i = 0; i < _n; i++) {
            result[i].reserve(group_size[i]);
        }
        for (int i = 0; i < _n; i++) {
            result[leader_buf[i]].push_back(i);
        }
        result.erase(
            std::remove_if(result.begin(), result.end(),
                           [&](const std::vector<int>& v) { return v.empty(); }),
            result.end());
        return result;
    }

  private:
    int _n;
    // root node: -1 * component size
    // otherwise: parent
    std::vector<int> parent_or_size;
};

// 連結グラフがループを持っているならtrueを返す
// それ以外はfalseを返す
bool has_loop(vector<vector<pair<int,double>>> edges){
    // 訪問済みのノードを記録する
    vector<int> visited(edges.size(),-1);

    // nodeとparentを持つスタック
    // DFSでループを検知する
    stack<pair<int,int>> stk;
    stk.push(pair<int,int>(0,-1));
    visited[0] = 1;
    while(!stk.empty()){
        auto [node,parent] = stk.top();
        stk.pop();
        for(auto [v,_] : edges[node]){
            if(v == parent) continue;
            if(visited[v] != -1) return true;
            visited[v] = 1;
            stk.push(pair<int,int>(v,node));
        }
    }
    return false;
}

int main() {
    // 辺を格納するための配列
    // 各ノードのインデックスに対して、隣接するノードとその重みを格納する
    vector<vector<pair<int,double>>> edges;

    // 駅のIDとインデックスを対応させるためのマップ
    // キーは駅のID、値はそのノードのインデックス
    // 駅のIDが1から始まると指定されていないため、用意する
    map<int,int> id_to_index;
    // インデックスと駅のIDを対応させるためのマップ
    // キーはノードのインデックス、値は駅のID
    map<int,int> index_to_id;
    
    // 入力を受け取る
    string line;
    int node_index = 0;
    while (getline(cin, line)) {
        stringstream ss(line);
        string s;
        int u, v;
        double w;

        getline(ss, s, ',');
        u = stoi(s);

        getline(ss, s, ',');
        v = stoi(s);

        getline(ss, s);
        w = stod(s);

        // ノードのインデックスがまだ存在しない場合、追加する
        if(id_to_index.count(u) == 0) {
            index_to_id[node_index] = u;
            id_to_index[u] = node_index++;
            edges.emplace_back(vector<pair<int,double>>());
        }
        if(id_to_index.count(v) == 0) {
            index_to_id[node_index] = v;
            id_to_index[v] = node_index++;
            edges.emplace_back(vector<pair<int,double>>());
        }

        // 辺を追加する
        edges[id_to_index[u]].emplace_back(pair<int,double>(id_to_index[v],w));
        edges[id_to_index[v]].emplace_back(pair<int,double>(id_to_index[u],w));
    }

    // union_findで連結成分を調べる
    // 連結でない成分同士は行き来できないので、連結なグラフごとに考えればよい
    dsu union_find(edges.size());
    for(int i = 0;i < edges.size();i++){
        for(int j = 0;j < edges[i].size();j++){
            union_find.merge(i,edges[i][j].first);
        }
    }

    // 連結なグラフごとの辺を格納するための配列
    vector<vector<vector<pair<int,double>>>> edges_groups;
    vector<vector<int>> groups = union_find.groups();
    for(int i = 0;i < groups.size();i++){
        edges_groups.emplace_back(vector<vector<pair<int,double>>>());

        map<int,int> global_to_local;
        for(int j = 0;j < groups[i].size();j++){
            global_to_local[groups[i][j]] = j; 
        }

        for(int j = 0;j < groups[i].size();j++){
            edges_groups[i].emplace_back(edges[groups[i][j]]);
            for(int k = 0;k < edges_groups[i][j].size();k++){
                edges_groups[i][j][k].first = global_to_local[edges_groups[i][j][k].first];
            }
        }
    }

    // 最長経路の距離と経路
    double max_distance = 0;
    vector<int> max_route;

    // 連結なグラフごとに最長経路を求める。
    for(int groups_idx = 0;groups_idx < edges_groups.size();groups_idx++){
        // ビットDPでO(2^N*N^2)で求める
        if(has_loop(edges_groups[groups_idx])){    
            // 最長な経路の情報を記録する
            // 距離
            double farthest_dist = 0.0;
            // 状態(訪問した点をbitで管理)
            long long farthest_state = -1;
            // 終点と始点
            // farthest_endはループなら終点の一個手前
            int farthest_end = -1,farthest_start = 1;
            // 始点と終点が一致しているかどうか？
            bool is_loop = false;
            
            // dp[S][v][p]はpからスタートして状態S(訪問済みノードをビットで管理したもの)で、vを最後に訪れたときの最長距離
            // pを持たないと1->2->3と2->1->3を同じものとして扱われ始点と終点が一致するケースに対応できなくなるのでpは必要
            vector<vector<vector<double>>> dp(1LL<<edges_groups[groups_idx].size(),vector<vector<double>>(edges_groups[groups_idx].size(),vector<double>(edges_groups[groups_idx].size(),-1.0)));
            // parent_prev[S][v][p]はpからスタートして状態S(訪問済みノードをビットで管理したもの)で、vを最後に訪れたときの一個前に訪問した点とスタート地点の点のペア
            vector<vector<vector<int>>> prev(1LL<<edges_groups[groups_idx].size(), vector<vector<int>>(edges_groups[groups_idx].size(), vector<int>(edges_groups[groups_idx].size(),-1)));

            // 初期化(各ノードから出発)
            for(long long i = 0;i < edges_groups[groups_idx].size();i++) {
                dp[1LL << i][i][i] = 0.0;
            }

            for(long long S = 0;S < 1LL<<edges_groups[groups_idx].size();S++){
                for(int v = 0;v < edges_groups[groups_idx].size();v++) {
                    for(int p = 0;p < edges_groups[groups_idx].size();p++){
                        // dp[S][v][p] < 0だとスタートがpで状態Sで点vから始まることができないのでcontinue
                        if(dp[S][v][p] < 0) continue;

                        // 点jに繋がる辺を調べる
                        for(int k = 0;k < edges_groups[groups_idx][v].size();k++){

                            // 始点に戻ってかつ最長ならその経路を記録
                            if(edges_groups[groups_idx][v][k].first == p && 
                            farthest_dist < dp[S][v][p]+edges_groups[groups_idx][v][k].second){
                                farthest_dist =  dp[S][v][p]+edges_groups[groups_idx][v][k].second;
                                farthest_state = S;
                                farthest_start = p;
                                farthest_end = v;
                                is_loop = true;
                                continue;
                            }

                            //すでに訪問済みな点はパス
                            if(S & (1LL<<edges_groups[groups_idx][v][k].first))continue;

                            // 経路が最長なら記録
                            if(dp[S|(1LL<<edges_groups[groups_idx][v][k].first)][edges_groups[groups_idx][v][k].first][p] < dp[S][v][p]+edges_groups[groups_idx][v][k].second){
                                dp[S|(1LL<<edges_groups[groups_idx][v][k].first)][edges_groups[groups_idx][v][k].first][p] = dp[S][v][p]+edges_groups[groups_idx][v][k].second;
                                prev[S|(1LL<<edges_groups[groups_idx][v][k].first)][edges_groups[groups_idx][v][k].first][p] = v;
                                if(farthest_dist < dp[S][v][p]+edges_groups[groups_idx][v][k].second){
                                    farthest_dist = dp[S][v][p]+edges_groups[groups_idx][v][k].second;
                                    farthest_state = S|(1LL<<edges_groups[groups_idx][v][k].first);
                                    farthest_start = p;
                                    farthest_end = edges_groups[groups_idx][v][k].first;
                                    is_loop = false;
                                }
                            }
                        }
                    }
                }
            }

            // 最長経路の距離の更新
            if(farthest_dist > max_distance){
                // 経路の復元
                vector<int> route;
                // 始点と終点が一致していたらを始点(終点)を追加
                if(is_loop)route.emplace_back(farthest_start);
                long long state = farthest_state, node = farthest_end;
                // それぞれの状態から一個前に訪問した点を順に探索することで経路を復元
                while (node != -1) {
                    route.push_back(node);
                    int prev_node = prev[state][node][farthest_start];
                    if (prev_node == -1) break;
                    state ^= (1LL << node);
                    node = prev_node;
                }
                
                // edges_groupsのインデックスからedgesのインデックスに復元
                for(int i = 0;i < route.size();i++){
                    route[i] = groups[groups_idx][route[i]];
                }

                max_distance = farthest_dist;
                max_route = route;
            }   
            
        }else{
            // ループがない(木構造)ならO(N)で最長経路を求めることができる。(Nは駅の個数)
            double farthest_dist = 0.0;
            vector<int> farthest_route;

            auto get_farthest = [&](int start) {
                // dfsでstartから一番遠い経路を求める
                // node,parent,distance,route
                stack<tuple<int,int,double,vector<int>>> stk;
                stk.push(tuple<int,int,double,vector<int>>(start,-1,0,vector<int>({start})));
                while(!stk.empty()){
                    auto [node, parent, sum_dist, route] = stk.top();
                    stk.pop();
                    for(auto [v,dist] : edges_groups[groups_idx][node]){
                        if(v == parent) continue;
                        vector<int> new_route = route;
                        new_route.emplace_back(v);
                        if(sum_dist + dist > farthest_dist){
                            //もっとも遠い距離、経路を記録
                            farthest_dist = sum_dist + dist;
                            farthest_route = new_route;
                        }
                        stk.push(tuple<int,int,double,vector<int>>(v,node,sum_dist + dist,new_route));
                    }
                }
            };

            // 0からもっとも遠い点を求める
            get_farthest(0);

            // 0からもっとも遠い点が直径(最長経路)の端点となる
            int diameter_start = farthest_route[farthest_route.size()-1];

            // 直径(最長経路)の端点からもっとも遠い点を見つけることで木の直径を求める
            get_farthest(diameter_start);

            // 木で始点と終点を同じにする場合は0->1->0のように往復する
            // 往復する場合では最長の辺を往復したものがもっとも長い
            for(int i = 0;i < edges_groups[groups_idx].size();i++){
                for(int j = 0;j < edges_groups[groups_idx][i].size();j++){
                    if(farthest_dist < 2 * edges_groups[groups_idx][i][j].second){
                        farthest_dist = 2 * edges_groups[groups_idx][i][j].second;
                        farthest_route = vector<int>{i,edges_groups[groups_idx][i][j].first,i};
                    }
                }
            }

            // この木の直径がすべての駅の経路の中で最長なら更新
            if(farthest_dist > max_distance){
                // edges_groupsのインデックスからedgesのインデックスに復元
                for(int i = 0;i < farthest_route.size();i++){
                    farthest_route[i] = groups[groups_idx][farthest_route[i]];
                }

                max_distance = farthest_dist;
                max_route = farthest_route;
            }
        }
    }
    
    // 経路を元のidに復元
    vector<int> max_route_id(max_route.size());
    for(int i = 0;i < max_route.size();i++){
        max_route_id[i] = index_to_id[max_route[i]];
    }

    // 見栄えのために始点のidが小さくなるようにする
    if(max_route[0] > max_route[max_route.size()])reverse(max_route_id.begin(),max_route_id.end());
    
    // 出力
    for(int i = 0;i < max_route_id.size();i++){
        cout << max_route_id[i] << endl;
    }

    return 0;
}