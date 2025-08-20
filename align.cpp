#include <algorithm>
#include <chrono>
#include <future>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <atomic>
#include <cctype>
#include <cstdlib>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::pair;

// check for valid sequences and exit if it is not
int valid_sequences(const string& s1, const string& s2);
void align_and_report(const string& s1, const string& s2,
                      bool score_only, bool one_alignment,
                      size_t max_alignments, int max_tasks);

// async traceback signature
std::vector<std::pair<std::string,std::string>>
traceback_all_optimal_mt(int i, int j,
                         const vector<vector<int>>& pointer_matrix,
                         const string& s1, const string& s2,
                         string buf1, string buf2,
                         size_t max_alignments,
                         std::atomic<int>& active_tasks,
                         int max_tasks);

// deterministic single traceback (greedy tie-break DIAG > UP > LEFT)
std::pair<string,string> deterministic_traceback(const vector<vector<int>>& P,
                                                const string& s1, const string& s2);

// count number of optimal alignments (ways) with cap
unsigned long long count_optimal_paths_capped(const vector<vector<int>>& P,
                                              int m, int n,
                                              unsigned long long cap);

// read FASTA or raw sequence from file: concatenate non-header lines, strip whitespace, uppercase
bool read_sequence_from_file(const string& path, string& out_seq, string& err_msg) {
  std::ifstream ifs(path);
  if (!ifs) {
    err_msg = "Cannot open file: " + path;
    return false;
  }
  std::string line;
  out_seq.clear();
  while (std::getline(ifs, line)) {
    if (!line.empty() && line[0] == '>') continue; // skip header
    for (char c : line) {
      if (!std::isspace(static_cast<unsigned char>(c))) {
        out_seq.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
      }
    }
  }
  if (out_seq.empty()) {
    err_msg = "No sequence data found in file: " + path;
    return false;
  }
  return true;
}

void print_usage(const char* prog) {
  cerr << "Usage: " << prog
       << " [--score-only] [--one-alignment] [--max-alignments N] [--max-tasks N]\n"
       << "           [--file1 <path> --file2 <path>] <sequence1> <sequence2>\n\n";
  cerr << "If --file1 and --file2 are provided, sequences are read from those files (FASTA or raw).\n";
  cerr << "If files are not provided, the last two positional arguments are treated as raw sequences.\n";
  cerr << "Defaults: max-alignments=10000, max-tasks=hardware_concurrency or 4\n";
}

int main(int argc, char* argv[]) {
  if (argc < 3) { print_usage(argv[0]); return 1; }

  bool score_only = false;
  bool one_alignment = false;
  size_t max_alignments = 10000;
  int max_tasks = 0; // if 0 -> auto (hw conc or 4)
  string file1, file2;

  // simple arg parsing
  int idx = 1;
  while (idx < argc) {
    string a = argv[idx];
    if (a == "--score-only") { score_only = true; ++idx; continue; }
    if (a == "--one-alignment") { one_alignment = true; ++idx; continue; }
    if (a == "--max-alignments") {
      if (idx + 1 >= argc) { cerr << "--max-alignments needs a number\n"; return 1; }
      max_alignments = std::stoull(argv[idx+1]);
      idx += 2; continue;
    }
    if (a == "--max-tasks") {
      if (idx + 1 >= argc) { cerr << "--max-tasks needs a number\n"; return 1; }
      max_tasks = std::stoi(argv[idx+1]);
      idx += 2; continue;
    }
    if (a == "--file1") {
      if (idx + 1 >= argc) { cerr << "--file1 needs a path\n"; return 1; }
      file1 = argv[idx+1];
      idx += 2; continue;
    }
    if (a == "--file2") {
      if (idx + 1 >= argc) { cerr << "--file2 needs a path\n"; return 1; }
      file2 = argv[idx+1];
      idx += 2; continue;
    }
    // otherwise break to collect positional args
    break;
  }

  // remaining positional args
  int remaining = argc - idx;
  string seq1_arg, seq2_arg;
  if (!file1.empty() || !file2.empty()) {
    // require both files if using file flags
    if (file1.empty() || file2.empty()) {
      cerr << "Both --file1 and --file2 must be provided when using file input.\n";
      return 1;
    }
    // read both files
    string err;
    if (!read_sequence_from_file(file1, seq1_arg, err)) { cerr << err << "\n"; return 1; }
    if (!read_sequence_from_file(file2, seq2_arg, err)) { cerr << err << "\n"; return 1; }
  } else {
    // expect two positional sequences
    if (remaining != 2) { print_usage(argv[0]); return 1; }
    seq1_arg = argv[idx];
    seq2_arg = argv[idx+1];
    // convert to uppercase and remove whitespace (defensive)
    string tmp;
    for (char c : seq1_arg) if (!std::isspace(static_cast<unsigned char>(c))) tmp.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    seq1_arg.swap(tmp);
    tmp.clear();
    for (char c : seq2_arg) if (!std::isspace(static_cast<unsigned char>(c))) tmp.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    seq2_arg.swap(tmp);
  }

  if (valid_sequences(seq1_arg, seq2_arg) != 0) return 1;

  // set default max_tasks if not provided
  if (max_tasks <= 0) {
    unsigned int hw = std::thread::hardware_concurrency();
    max_tasks = (hw > 0 ? static_cast<int>(hw) : 4);
  }

  align_and_report(seq1_arg, seq2_arg, score_only, one_alignment, max_alignments, max_tasks);
  return 0;
}

int valid_sequences(const string& seq1, const string& seq2) {
  for (const char n : seq1) {
    if (n != 'A' && n != 'G' && n != 'C' && n != 'T') {
      cerr << "Invalid sequence (seq1 contains non-ACGT characters)\n";
      return 1;
    }
  }
  for (const char n : seq2) {
    if (n != 'A' && n != 'G' && n != 'C' && n != 'T') {
      cerr << "Invalid sequence (seq2 contains non-ACGT characters)\n";
      return 1;
    }
  }
  if (seq1.empty() || seq2.empty()) {
    cerr << "There is an empty sequence\n";
    return 1;
  }
  return 0;
}

// Main function: build matrices, optionally count ways, optionally enumerate
void align_and_report(const string& s1, const string& s2,
                      bool score_only, bool one_alignment,
                      size_t max_alignments, int max_tasks) {
  const int match = 1;
  const int mismatch = 1;
  const int gap = 1;
  const int DIAG = 1;
  const int UP = 2;
  const int LEFT = 4;

  int m = (int)s1.length();
  int n = (int)s2.length();

  // start timing DP
  using clk = std::chrono::high_resolution_clock;
  auto t0 = clk::now();

  vector<vector<int>> S(m + 1, vector<int>(n + 1, -1));
  vector<vector<int>> P(m + 1, vector<int>(n + 1, 0));

  for (int i = 0; i <= m; ++i) S[i][0] = -i;
  for (int j = 0; j <= n; ++j) S[0][j] = -j;
  for (int i = 1; i <= m; ++i) P[i][0] = UP;
  for (int j = 1; j <= n; ++j) P[0][j] = LEFT;
  P[0][0] = 0;

  for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; ++j) {
      int up   = S[i-1][j] - gap;
      int lf   = S[i][j-1] - gap;
      int diag = S[i-1][j-1] + ((s1[i-1] == s2[j-1]) ? match : -mismatch);
      int best = std::max(std::max(diag, up), lf);
      S[i][j] = best;
      int mask = 0;
      if (best == diag) mask |= DIAG;
      if (best == up)   mask |= UP;
      if (best == lf)   mask |= LEFT;
      P[i][j] = mask;
    }
  }

  auto t1 = clk::now();
  long long dp_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

  // Print final score
  cout << S[m][n] << endl;

  // If user only wanted the score:
  if (score_only) {
    cerr << "DP fill time: " << dp_ms << " ms\n";
    return;
  }

  // Count number of optimal alignments (capped at max_alignments)
  unsigned long long cap = (max_alignments > 0 ? max_alignments : 10000ULL);
  unsigned long long ways = count_optimal_paths_capped(P, m, n, cap);

  if (ways >= cap) {
    cerr << "Number of optimal alignments >= " << cap << " (capped).\n";
    if (!one_alignment) {
      cerr << "Enumeration skipped to avoid explosion. Use --one-alignment to get one alignment, or increase --max-alignments.\n";
      cerr << "DP fill time: " << dp_ms << " ms\n";
      return;
    } else {
      // proceed to deterministic single traceback
      auto single = deterministic_traceback(P, s1, s2);
      cout << single.first << endl << single.second << endl;
      cerr << "Number of optimal alignments >= " << cap << " (capped). Returned one alignment (deterministic).\n";
      cerr << "DP fill time: " << dp_ms << " ms\n";
      return;
    }
  }

  // ways < cap -> safe to enumerate (up to max_alignments)
  cerr << "Number of optimal alignments = " << ways << " (will enumerate up to " << max_alignments << ")\n";

  // Prepare async traceback
  std::atomic<int> active_tasks(0);
  auto ttrace_start = clk::now();
  vector<std::pair<std::string,std::string>> allAlignments =
    traceback_all_optimal_mt(m, n, P, s1, s2, string(""), string(""),
                             max_alignments, active_tasks, max_tasks);
  auto ttrace_end = clk::now();
  long long trace_ms = std::chrono::duration_cast<std::chrono::milliseconds>(ttrace_end - ttrace_start).count();

  // print alignments (main thread)
  for (const auto &p : allAlignments) {
    cout << p.first << endl;
    cout << p.second << endl;
    cout << endl;
  }

  if (allAlignments.size() >= max_alignments) {
    cerr << "Returned alignments were truncated at max_alignments = " << max_alignments << "\n";
  }
  cerr << "DP fill time: " << dp_ms << " ms\n";
  cerr << "Traceback/enumeration time: " << trace_ms << " ms\n";
}

// deterministic greedy single traceback (DIAG > UP > LEFT)
std::pair<string,string> deterministic_traceback(const vector<vector<int>>& P,
                                                const string& s1, const string& s2) {
  int m = (int)s1.size();
  int n = (int)s2.size();
  string a1, a2;
  int i = m, j = n;
  while (i > 0 || j > 0) {
    int mask = P[i][j];
    if ((mask & 1) && i > 0 && j > 0) {
      a1.push_back(s1[i-1]);
      a2.push_back(s2[j-1]);
      --i; --j;
    } else if ((mask & 2) && i > 0) {
      a1.push_back(s1[i-1]);
      a2.push_back('-');
      --i;
    } else if ((mask & 4) && j > 0) {
      a1.push_back('-');
      a2.push_back(s2[j-1]);
      --j;
    } else {
      break;
    }
  }
  std::reverse(a1.begin(), a1.end());
  std::reverse(a2.begin(), a2.end());
  return {a1, a2};
}

// Count how many optimal paths exist using P; cap counts at 'cap'
unsigned long long count_optimal_paths_capped(const vector<vector<int>>& P,
                                              int m, int n,
                                              unsigned long long cap) {
  vector<vector<unsigned long long>> ways(m+1, vector<unsigned long long>(n+1, 0));
  ways[0][0] = 1;
  for (int i = 0; i <= m; ++i) {
    for (int j = 0; j <= n; ++j) {
      if (i==0 && j==0) continue;
      unsigned long long w = 0;
      int mask = P[i][j];
      if ((mask & 1) && i-1 >= 0 && j-1 >= 0) { w += ways[i-1][j-1]; if (w >= cap) { ways[i][j] = cap; continue; } }
      if ((mask & 2) && i-1 >= 0) { w += ways[i-1][j]; if (w >= cap) { ways[i][j] = cap; continue; } }
      if ((mask & 4) && j-1 >= 0) { w += ways[i][j-1]; if (w >= cap) { ways[i][j] = cap; continue; } }
      ways[i][j] = (w >= cap ? cap : w);
    }
  }
  return ways[m][n];
}

// async-enabled traceback (returns vector of alignments). Throttles spawn using active_tasks and max_tasks.
std::vector<std::pair<std::string,std::string>>
traceback_all_optimal_mt(int i, int j,
                         const vector<vector<int>>& pointer_matrix,
                         const string& s1, const string& s2,
                         string buf1, string buf2,
                         size_t max_alignments,
                         std::atomic<int>& active_tasks,
                         int max_tasks) {
  const int DIAG = 1, UP = 2, LEFT = 4;
  using ResT = std::vector<std::pair<std::string,std::string>>;
  ResT results;

  if (i == 0 && j == 0) {
    std::reverse(buf1.begin(), buf1.end());
    std::reverse(buf2.begin(), buf2.end());
    results.emplace_back(std::move(buf1), std::move(buf2));
    return results;
  }

  int mask = pointer_matrix[i][j];
  std::vector<std::future<ResT>> futures;

  auto spawn_or_run = [&](auto work_lambda) {
    int prev = active_tasks.fetch_add(1);
    if (prev < max_tasks) {
      auto fut = std::async(std::launch::async, [work_lambda, &active_tasks]() -> ResT {
        ResT r = work_lambda();
        active_tasks.fetch_sub(1);
        return r;
      });
      futures.push_back(std::move(fut));
    } else {
      active_tasks.fetch_sub(1);
      ResT r = work_lambda();
      for (auto &p : r) {
        if (results.size() >= max_alignments) return;
        results.push_back(std::move(p));
      }
    }
  };

  bool first_branch_done = false;

  // DIAG
  if (mask & DIAG) {
    auto work = [=, &pointer_matrix, &s1, &s2, &active_tasks]() -> ResT {
      string nb1 = buf1; nb1.push_back(s1[i-1]);
      string nb2 = buf2; nb2.push_back(s2[j-1]);
      return traceback_all_optimal_mt(i-1, j-1, pointer_matrix, s1, s2,
                                      std::move(nb1), std::move(nb2),
                                      max_alignments, active_tasks, max_tasks);
    };
    if (!first_branch_done) {
      ResT r = work();
      for (auto &p : r) {
        if (results.size() >= max_alignments) return results;
        results.push_back(std::move(p));
      }
      first_branch_done = true;
    } else {
      spawn_or_run(work);
    }
  }

  // UP
  if (mask & UP) {
    auto work = [=, &pointer_matrix, &s1, &s2, &active_tasks]() -> ResT {
      string nb1 = buf1; nb1.push_back(s1[i-1]);
      string nb2 = buf2; nb2.push_back('-');
      return traceback_all_optimal_mt(i-1, j, pointer_matrix, s1, s2,
                                      std::move(nb1), std::move(nb2),
                                      max_alignments, active_tasks, max_tasks);
    };
    if (!first_branch_done) {
      ResT r = work();
      for (auto &p : r) {
        if (results.size() >= max_alignments) return results;
        results.push_back(std::move(p));
      }
      first_branch_done = true;
    } else {
      spawn_or_run(work);
    }
  }

  // LEFT
  if (mask & LEFT) {
    auto work = [=, &pointer_matrix, &s1, &s2, &active_tasks]() -> ResT {
      string nb1 = buf1; nb1.push_back('-');
      string nb2 = buf2; nb2.push_back(s2[j-1]);
      return traceback_all_optimal_mt(i, j-1, pointer_matrix, s1, s2,
                                      std::move(nb1), std::move(nb2),
                                      max_alignments, active_tasks, max_tasks);
    };
    if (!first_branch_done) {
      ResT r = work();
      for (auto &p : r) {
        if (results.size() >= max_alignments) return results;
        results.push_back(std::move(p));
      }
      first_branch_done = true;
    } else {
      spawn_or_run(work);
    }
  }

  // collect futures
  for (auto &f : futures) {
    ResT branch_res = f.get();
    for (auto &p : branch_res) {
      if (results.size() >= max_alignments) return results;
      results.push_back(std::move(p));
    }
  }

  return results;
}
