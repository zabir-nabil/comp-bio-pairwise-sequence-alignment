#include <bits/stdc++.h>

#ifdef __linux__
#   include <sys/resource.h>   // getrusage
static long peak_rss_MB() {
    struct rusage ru{};
    getrusage(RUSAGE_SELF, &ru);
    // ru.ru_maxrss is kilobytes on Linux → convert to MB
    return ru.ru_maxrss / 1024;
}
#elif _WIN32
#   include <windows.h>
#   include <psapi.h>
static long peak_rss_MB() {
    PROCESS_MEMORY_COUNTERS info{};
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return static_cast<long>(info.PeakWorkingSetSize / (1024*1024));
}
#else
static long peak_rss_MB() { return -1; } // unsupported OS
#endif


using namespace std;

/* ------------------------------------------------------------- */
/*  i store the scoring matrix and the DNA DNA_GAP_ALPHABET here */
/* ------------------------------------------------------------- */
static int  score_matrix[5][5];              
static const string DNA_GAP_ALPHABET = "ACGT-";

/* ------------------------------------------------------------- */
/*  i read an entire text file and squeeze it into one long line  */
/*  every whitespace character (\n, \r, spaces, tabs) disappears */
/*  so the caller gets a clean sequence string in return          */
/* ------------------------------------------------------------- */
static string load_sequence(const string& path)
{
    ifstream in(path);
    if (!in) {
        throw runtime_error("couldn't open " + path);
    }

    string line, sequence;
    while (getline(in, line)) {
        /* delete all leading+trailing whitespace */
        line.erase(remove_if(line.begin(), line.end(),
                             [](unsigned char c){ return isspace(c); }),
                   line.end());
        sequence += line;   // append the cleaned‑up chunk
    }
    return sequence;
}


/* -------------------------------------------------------------- */
/*  Rotate 2-D matrix 180°                                        */
/* -------------------------------------------------------------- */
static void rotate_matrix(vector<vector<int>>& M)
{
    const int rows = (int)M.size();
    const int cols = (int)M[0].size();

    /* reverse each row */
    for (auto& row : M) std::reverse(row.begin(), row.end());
    /* reverse order of rows */
    std::reverse(M.begin(), M.end());
}

/* ------------------------------------------------------------- */
/*  i compute a 2‑D dynamic‑programming table for two slices      */
/*  of the sequences – this is the pairwise sub‑problem that the  */
/*  big 3‑D alignment reuses again and again.                     */
/*                                                               */
/*  the routine is a thin wrapper around Needleman‑Wunsch:        */
/*  - rows   == slice of string x (xBeg…xEnd)                     */
/*  - cols   == slice of string y (yBeg…yEnd)                     */
/*  - gap penalties are linear (always 2 * score_matrix[4][*])    */
/*  i keep it separate so i can test it in isolation.             */
/* ------------------------------------------------------------- */
static vector<vector<int>>
score2d_region(const string& x,const string& y,
               int xBeg,int yBeg,int xEnd,int yEnd)
{
    const int rows = xEnd - xBeg;
    const int cols = yEnd - yBeg;

    /* dp[r][c] = best score aligning first r symbols of x‑slice   */
    /*           with first c symbols of y‑slice                   */
    vector<vector<int>> dp(rows+1, vector<int>(cols+1, 0));

    // first column : x‑slice against all gaps in y
    for (int r = 1; r <= rows; ++r) {
        int idx = DNA_GAP_ALPHABET.find( x[xBeg + r - 1] );
        dp[r][0] = dp[r-1][0] + 2 * score_matrix[4][idx];
    }
    // first row : y‑slice against all gaps in x
    for (int c = 1; c <= cols; ++c) {
        int idx = DNA_GAP_ALPHABET.find( y[yBeg + c - 1] );
        dp[0][c] = dp[0][c-1] + 2 * score_matrix[4][idx];
    }

    /* interior cells – classic three‑way max                       */
    for (int r = 1; r <= rows; ++r) {
        int idxX = DNA_GAP_ALPHABET.find( x[xBeg + r - 1] );
        for (int c = 1; c <= cols; ++c) {
            int idxY = DNA_GAP_ALPHABET.find( y[yBeg + c - 1] );

            int diag = dp[r-1][c-1] +                   // match / mismatch
                       score_matrix[idxX][idxY] +      // pair score
                       score_matrix[4][idxX] +         // gap cost contributions
                       score_matrix[4][idxY];

            int up   = dp[r-1][c]   + 2 * score_matrix[4][idxX]; // gap in y
            int left = dp[r][c-1]   + 2 * score_matrix[4][idxY]; // gap in x

            dp[r][c] = std::max({diag, up, left});
        }
    }

    return dp;
}


/* ------------------------------------------------------------------ */
/*  frontier3d_plane                                                   */
/*                                                                    */
/*  builds one forward "slab" of the 3‑way alignment cube.            */
/*                                                                    */
/*  imagine a 3‑D DP box indexed by (i, j, k) where                    */
/*    i walks along sequence x                                         */
/*    j walks along sequence y                                         */
/*    k walks along sequence z                                         */
/*                                                                    */
/*  we only need a sub‑box:                                            */
/*    i in [xb, xe)   – start (xb) inclusive, end (xe) exclusive       */
/*    j in [yb, ye)                                                   */
/*    k in (zb, ze]   – we extend layer‑by‑layer from right after zb   */
/*                       up to and INCLUDING ze.                       */
/*  (so dz = ze ‑ zb is the number of z‑layers we will add)            */
/*                                                                    */
/*  the function returns a 2‑D matrix F of size (dx+1) × (dy+1) where   */
/*    F[i][j] holds the best score ending at                           */
/*              (i + xb,  j + yb,  current k)                          */
/*  only O(dx·dy) memory is used by re‑using the buffer for each k.    */
/*                                                                    */
/*  parameters                                                         */
/*    x, y, z : the full sequences                                     */
/*    xb,xe   : slice of x (begin, end)                                */
/*    yb,ye   : slice of y                                             */
/*    zb,ze   : z‑range we sweep through                               */
/* ------------------------------------------------------------------ */
static vector<vector<int>>
frontier3d_plane(const string& x,const string& y,const string& z,
                 int xb,int yb,int zb,
                 int xe,int ye,int ze)
{
    /* pre‑compute pairwise tables for the whole sub‑box */
    auto xy = score2d_region(x,y, xb,yb, xe,ye);
    auto zx = score2d_region(z,x, zb,xb, ze,xe);
    auto zy = score2d_region(z,y, zb,yb, ze,ye);

    int dx = xe - xb;      // length of x slice
    int dy = ye - yb;      // length of y slice
    int dz = ze - zb;      // how many z layers we will add

    vector<vector<int>> F = xy;   // workspace for the current z layer

    for (int k=1; k<=dz; ++k) {
        /* 1. move previous layer into xy (we need it for transitions) */
        for (int i=0;i<=dx;++i)
            for (int j=0;j<=dy;++j)
                xy[i][j] = F[i][j];

        /* 2. clear F so we can fill the new layer */
        for (auto& row : F) fill(row.begin(), row.end(), 0);

        /* 3. set borders – aligning x or y fully with gaps */
        for (int i=0;i<=dx;++i) F[i][0] = zx[k][i];
        for (int j=0;j<=dy;++j) F[0][j] = zy[k][j];

        /* 4. interior cells: evaluate 7 possibilities */
        for (int i=1;i<=dx;++i)
            for (int j=1;j<=dy;++j) {
                int a = DNA_GAP_ALPHABET.find( x[xb+i-1] );
                int b = DNA_GAP_ALPHABET.find( y[yb+j-1] );
                int c = DNA_GAP_ALPHABET.find( z[zb+k-1] );

                int p00 = xy[i-1][j-1]; // all from previous layer
                int p10 = xy[i-1][j];
                int p01 = xy[i][j-1];
                int p11 = xy[i][j];

                int q11 = F[i-1][j-1]; // already in this layer
                int q01 = F[i][j-1];
                int q10 = F[i-1][j];

                int gap_x   = q10 + 2*score_matrix[4][a];
                int gap_y   = q01 + 2*score_matrix[4][b];
                int gap_z   = p11 + 2*score_matrix[4][c];
                int pair_xy = q11 + score_matrix[4][a] + score_matrix[4][b] + score_matrix[a][b];
                int pair_xz = p10 + score_matrix[4][a] + score_matrix[4][c] + score_matrix[a][c];
                int pair_yz = p01 + score_matrix[4][b] + score_matrix[4][c] + score_matrix[b][c];
                int trio    = p00 + score_matrix[a][b] + score_matrix[a][c] + score_matrix[b][c];

                F[i][j] = max({gap_x, gap_y, gap_z, pair_xy, pair_xz, pair_yz, trio});
            }
    }
    return F;
}


/* dna 3‑way alignment */

/* ------------------------------------------------------------- */
/*  reverse_frontier_plane                                        */
/*  i run the forward frontier code on the reversed sequences     */
/*  then rotate the plane 180° so that indices line up with the    */
/*  original (non‑reversed) coordinate system.                    */
/* ------------------------------------------------------------- */
static vector<vector<int>>
reverse_frontier_plane(const string& rx,const string& ry,const string& rz,
                       int xb,int yb,int zb,
                       int xe,int ye,int ze)
{
    auto B = frontier3d_plane(rx, ry, rz, xb, yb, zb, xe, ye, ze);
    rotate_matrix(B); // flip so that (0,0) remains (0,0) in forward view
    return B;
}

/* ------------------------------------------------------------------
   cube_dp_align
   brute-force dynamic programme for a small alignment cube
   — i call the three local slices a, b, c
   — ma = a.size(), nb = b.size(), lc = c.size()
   — i keep one 3-D table ‘dp’ for scores and one ‘step’ for back-pointers
   ------------------------------------------------------------------ */

struct MiniCubeOut {
    int score;              // best score inside this cube
    string ra, rb, rc;      // aligned substrings (same length)
};

static MiniCubeOut cube_dp_align(const string& a,
                                 const string& b,
                                 const string& c)
{
    const int ma = a.size(), nb = b.size(), lc = c.size();

    /* dp[x][y][z]  = best score up to prefix (x,y,z)
       step[x][y][z] remembers *which* of the 7 moves was chosen        */
    vector dp(ma + 1,
              vector(nb + 1,
                     vector<int>(lc + 1, 0)));
    vector step(ma + 1,
                vector(nb + 1,
                       vector<int>(lc + 1, -1)));

    auto idx = [&](char ch){ return DNA_GAP_ALPHABET.find(ch); };

    /* --- fill edges where two axes are zero ------------------------ */
    for (int x = 1; x <= ma; ++x) {                // a-only axis
        dp[x][0][0] = dp[x-1][0][0] + 2*score_matrix[4][ idx(a[x-1]) ];
        step[x][0][0] = 0;                         // came from (x-1,0,0)
    }
    for (int y = 1; y <= nb; ++y) {                // b-only axis
        dp[0][y][0] = dp[0][y-1][0] + 2*score_matrix[4][ idx(b[y-1]) ];
        step[0][y][0] = 1;                         // came from (0,y-1,0)
    }
    for (int z = 1; z <= lc; ++z) {                // c-only axis
        dp[0][0][z] = dp[0][0][z-1] + 2*score_matrix[4][ idx(c[z-1]) ];
        step[0][0][z] = 2;                         // came from (0,0,z-1)
    }

    /* --- fill faces where exactly one axis is zero ----------------- */
    // xy-plane (z == 0)
    for (int x = 1; x <= ma; ++x)
        for (int y = 1; y <= nb; ++y) {
            int ia = idx(a[x-1]), ib = idx(b[y-1]);
            int diag = dp[x-1][y-1][0] + score_matrix[ia][ib] +
                        score_matrix[4][ia] + score_matrix[4][ib];
            int gapA = dp[x-1][y  ][0] + 2*score_matrix[4][ia];
            int gapB = dp[x  ][y-1][0] + 2*score_matrix[4][ib];
            dp[x][y][0] = max({diag, gapA, gapB});
            step[x][y][0] = (dp[x][y][0] == diag ? 3 :
                             dp[x][y][0] == gapA ? 0 : 1);
        }

    // xz-plane (y == 0)
    for (int x = 1; x <= ma; ++x)
        for (int z = 1; z <= lc; ++z) {
            int ia = idx(a[x-1]), ic = idx(c[z-1]);
            int diag = dp[x-1][0][z-1] + score_matrix[ia][ic] +
                        score_matrix[4][ia] + score_matrix[4][ic];
            int gapA = dp[x-1][0][z  ] + 2*score_matrix[4][ia];
            int gapC = dp[x  ][0][z-1] + 2*score_matrix[4][ic];
            dp[x][0][z] = max({diag, gapA, gapC});
            step[x][0][z] = (dp[x][0][z] == diag ? 4 :
                             dp[x][0][z] == gapA ? 0 : 2);
        }

    // yz-plane (x == 0)
    for (int y = 1; y <= nb; ++y)
        for (int z = 1; z <= lc; ++z) {
            int ib = idx(b[y-1]), ic = idx(c[z-1]);
            int diag = dp[0][y-1][z-1] + score_matrix[ib][ic] +
                        score_matrix[4][ib] + score_matrix[4][ic];
            int gapB = dp[0][y-1][z  ] + 2*score_matrix[4][ib];
            int gapC = dp[0][y  ][z-1] + 2*score_matrix[4][ic];
            dp[0][y][z] = max({diag, gapB, gapC});
            step[0][y][z] = (dp[0][y][z] == diag ? 5 :
                             dp[0][y][z] == gapB ? 1 : 2);
        }

    /* --- interior: genuine 3-D cells ------------------------------- */
    for (int z = 1; z <= lc; ++z)
        for (int x = 1; x <= ma; ++x)
            for (int y = 1; y <= nb; ++y) {

                int ia = idx(a[x-1]), ib = idx(b[y-1]), ic = idx(c[z-1]);

                int gapA = dp[x-1][y  ][z  ] + 2*score_matrix[4][ia];
                int gapB = dp[x  ][y-1][z  ] + 2*score_matrix[4][ib];
                int gapC = dp[x  ][y  ][z-1] + 2*score_matrix[4][ic];

                int ab   = dp[x-1][y-1][z  ] + score_matrix[4][ia] +
                            score_matrix[4][ib] + score_matrix[ia][ib];
                int ac   = dp[x-1][y  ][z-1] + score_matrix[4][ia] +
                            score_matrix[4][ic] + score_matrix[ia][ic];
                int bc   = dp[x  ][y-1][z-1] + score_matrix[4][ib] +
                            score_matrix[4][ic] + score_matrix[ib][ic];

                int abc  = dp[x-1][y-1][z-1] + score_matrix[ia][ib] +
                            score_matrix[ia][ic] + score_matrix[ib][ic];

                dp[x][y][z] = max({gapA, gapB, gapC, ab, ac, bc, abc});

                /* record which transition we took (0-6) */
                if      (dp[x][y][z] == gapA) step[x][y][z] = 0;
                else if (dp[x][y][z] == gapB) step[x][y][z] = 1;
                else if (dp[x][y][z] == gapC) step[x][y][z] = 2;
                else if (dp[x][y][z] == ab  ) step[x][y][z] = 3;
                else if (dp[x][y][z] == ac  ) step[x][y][z] = 4;
                else if (dp[x][y][z] == bc  ) step[x][y][z] = 5;
                else                          step[x][y][z] = 6;
            }

    /* --- traceback: build the aligned strings in reverse ----------- */
    string ra, rb, rc;
    int x = ma, y = nb, z = lc;

    while (x || y || z) {
        int s = step[x][y][z];          // remember the chosen move
        if (s == 0) {                   // gap in b and c
            ra += a[x-1]; rb += '-'; rc += '-'; --x;
        }
        else if (s == 1) {              // gap in a and c
            ra += '-'; rb += b[y-1]; rc += '-'; --y;
        }
        else if (s == 2) {              // gap in a and b
            ra += '-'; rb += '-'; rc += c[z-1]; --z;
        }
        else if (s == 3) {              // align a+b ; gap c
            ra += a[x-1]; rb += b[y-1]; rc += '-'; --x; --y;
        }
        else if (s == 4) {              // align a+c ; gap b
            ra += a[x-1]; rb += '-'; rc += c[z-1]; --x; --z;
        }
        else if (s == 5) {              // align b+c ; gap a
            ra += '-'; rb += b[y-1]; rc += c[z-1]; --y; --z;
        }
        else if (s == 6) {              // align a+b+c
            ra += a[x-1]; rb += b[y-1]; rc += c[z-1]; --x; --y; --z;
        }
        else {                          // sanity check
            throw logic_error("bad back-pointer");
        }
    }

    reverse(ra.begin(), ra.end());
    reverse(rb.begin(), rb.end());
    reverse(rc.begin(), rc.end());

    return { dp[ma][nb][lc], ra, rb, rc };
}


/* ------------------------------------------------------------------
   divide_conquer_align
   i split the 3-D problem cube in half along the c-sequence (z-axis).
   – if any edge of the cube is ≤1, i fall back to the cubic DP.
   – otherwise i score the “front” half, the “back” half (on reversed
     strings), find the best (x*,y*) cut, then recurse on both halves.
   the stitched fragments accumulate in the ‘out’ struct.
   ------------------------------------------------------------------ */
struct AlignOut {
    string sa, sb, sc;          // growing global alignment
    vector<int> frag_scores;    // score of each stitched sub-block
};

static void divide_conquer_align(            // divide and conquer
    const string& a,  const string& b,  const string& c,
    const string& ra, const string& rb, const string& rc,   // reversed
    int aBeg, int bBeg, int cBeg,                          // begin idx (inclusive)
    int aEnd, int bEnd, int cEnd,                          // end idx (exclusive)
    AlignOut& out)
{
    /* base-case: if any side ≤1 char, brute-force the little cube */
    if (aEnd - aBeg <= 1 || bEnd - bBeg <= 1 || cEnd - cBeg <= 1) {
        MiniCubeOut block = cube_dp_align(
            a.substr(aBeg, aEnd - aBeg),
            b.substr(bBeg, bEnd - bBeg),
            c.substr(cBeg, cEnd - cBeg));
        out.frag_scores.push_back(block.score);
        out.sa += block.ra; out.sb += block.rb; out.sc += block.rc;
        return;
    }

    /* split along the middle of the c-sequence ------------------- */
    int midZ = (cBeg + cEnd) / 2;

    /* scores from origin → (x,y,midZ) */
    auto front = frontier3d_plane(a, b, c,
                                  aBeg, bBeg, cBeg,
                                  aEnd, bEnd, midZ);

    /* map slice to reversed coordinates, then score
       (x,y,midZ) → terminus on reversed cube */
    int raBeg = int(a.size()) - aEnd, raEnd = int(a.size()) - aBeg;
    int rbBeg = int(b.size()) - bEnd, rbEnd = int(b.size()) - bBeg;
    int rcBeg = int(c.size()) - cEnd, rcEnd = int(c.size()) - midZ;

    auto back = reverse_frontier_plane(ra, rb, rc,
                                       raBeg, rbBeg, rcBeg,
                                       raEnd, rbEnd, rcEnd);

    /* choose cut point (x*,y*) that maximises front+back ---------- */
    int bestScore = INT_MIN, cutX = 0, cutY = 0;
    for (int x = 0; x < (int)front.size(); ++x)
        for (int y = 0; y < (int)front[x].size(); ++y) {
            int candidate = front[x][y] + back[x][y];
            if (candidate > bestScore) { bestScore = candidate; cutX = x; cutY = y; }
        }

    /* recurse on the two sub-boxes --------------------------------*/
    divide_conquer_align(a, b, c, ra, rb, rc,
                         aBeg,           bBeg,           cBeg,
                         aBeg + cutX,    bBeg + cutY,    midZ,
                         out);

    divide_conquer_align(a, b, c, ra, rb, rc,
                         aBeg + cutX,    bBeg + cutY,    midZ,
                         aEnd,           bEnd,           cEnd,
                         out);
}

/* -------------------------------------------------------------- */
/*  Initialise BLAST score matrix                                 */
/* -------------------------------------------------------------- */
static void setup_blast_score_matrix()
{
    score_matrix[0][0]=5;  score_matrix[0][1]=-4; score_matrix[0][2]=-4; score_matrix[0][3]=-4; score_matrix[0][4]=-8;
    score_matrix[1][0]=-4; score_matrix[1][1]=5;  score_matrix[1][2]=-4; score_matrix[1][3]=-4; score_matrix[1][4]=-8;
    score_matrix[2][0]=-4; score_matrix[2][1]=-4; score_matrix[2][2]=5;  score_matrix[2][3]=-4; score_matrix[2][4]=-8;
    score_matrix[3][0]=-4; score_matrix[3][1]=-4; score_matrix[3][2]=-4; score_matrix[3][3]=5;  score_matrix[3][4]=-8;
    score_matrix[4][0]=-8; score_matrix[4][1]=-8; score_matrix[4][2]=-8; score_matrix[4][3]=-8; score_matrix[4][4]=0;
}

/* -------------------------------------------------------------- */
/*                   =====     MAIN     =====                     */
/* -------------------------------------------------------------- */
int main()
{
    try{
        setup_blast_score_matrix();

        const string root = "./datasets/";
        
        /* for testing */
        // const string seq1 = root + "t1_s1.txt";
        // const string seq2 = root + "t1_s2.txt";
        // const string seq3 = root + "t1_s3.txt";
        
        
        const string seq1 = root + "NM_000492.txt";
        const string seq2 = root + "NM_021050.txt";
        const string seq3 = root + "NM_031506.txt";



        string s1 = load_sequence(seq1);
        string s2 = load_sequence(seq2);
        string s3 = load_sequence(seq3);

        cout<<"seq1 ("<< seq1 <<") length: "<<s1.size()<<'\n';
        cout<<"seq2 ("<< seq2 <<") length: "<<s2.size()<<'\n';
        cout<<"seq3 ("<< seq3 <<") length: "<<s3.size()<<'\n';

        string rs1(s1.rbegin(),s1.rend());
        string rs2(s2.rbegin(),s2.rend());
        string rs3(s3.rbegin(),s3.rend());

        AlignOut aln;
        auto t0 = chrono::high_resolution_clock::now();
        divide_conquer_align(s1,s2,s3, rs1,rs2,rs3,
                    0,0,0, (int)s1.size(), (int)s2.size(), (int)s3.size(),
                    aln);
        auto t1 = chrono::high_resolution_clock::now();
        cerr<<"Total running time (seconds): "
            << chrono::duration<double>(t1-t0).count() <<'\n' << "Peak RSS ≈ " << peak_rss_MB()
     << " MB\n"; // resident set size (RSS) is the portion of memory (measured in kilobytes) occupied by a process that is held in main memory (RAM)

        /* compute score again to cross-check */
        int score=0, matches=0;
        for(size_t i=0;i<aln.sa.size();++i){
            int ia=DNA_GAP_ALPHABET.find(aln.sa[i]);
            int ib=DNA_GAP_ALPHABET.find(aln.sb[i]);
            int ic=DNA_GAP_ALPHABET.find(aln.sc[i]);
            score += score_matrix[ia][ib] + score_matrix[ia][ic] + score_matrix[ib][ic];
            if (aln.sa[i]!='-' && aln.sa[i]==aln.sb[i] && aln.sa[i]==aln.sc[i]) ++matches;
        }

        cout<<"Alignment length of the sequences : "<<aln.sa.size()<<"\n";
        cout<<"# of exact matches  (identical XYZ) : "<<matches<<"\n";
        cout<<"Total score : "<<score<<"\n";
        // If want to print the raw alignment -> uncomment:
        cout<<aln.sa<<"\n"<<aln.sb<<"\n"<<aln.sc<<"\n";
        //
    }
    catch(const exception& ex){
        cerr<<"Error: "<<ex.what()<<"\n";
        return 1;
    }
    return 0;
}
