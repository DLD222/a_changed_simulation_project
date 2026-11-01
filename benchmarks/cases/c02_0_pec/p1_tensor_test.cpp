/* Standalone check of SubpixelTensor.H against p1_tensor_reference.py.
 *
 * Reads p1_cases.csv (case,kind,e1[9],e2[9],n[3],f), calls chi_inv_anisotropic
 * for every case (same dispatch as the Python reference), and writes p1_cpp.csv
 * in the format expected by `p1_tensor_reference.py --compare`.
 *
 * Build (does not touch the production target):
 *   g++ -O2 -std=c++17 -DAMREX_USE_ASSERTIONS \
 *       -I Source/FieldSolver/FiniteDifferenceSolver/MacroscopicProperties \
 *       -I build_2d/_deps/fetchedamrex-src/Src/Base \
 *       -I build_2d/_deps/fetchedamrex-src/Src/Base/Parser \
 *       p1_tensor_test.cpp -o /tmp/p1_tensor_test
 */
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "SubpixelTensor.H"

namespace {

std::vector<std::string> split(const std::string& line)
{
    std::vector<std::string> out;
    std::string cell;
    std::istringstream ss(line);
    while (std::getline(ss, cell, ',')) { out.push_back(cell); }
    return out;
}

// The CSV carries a full 3x3; the reference symmetrises it, and so do we.
warpx::subpixel::Mat3 sym_from_row(const std::vector<std::string>& v, int offset)
{
    amrex::Real a[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            a[i][j] = std::stod(v[offset + 3 * i + j]);
        }
    }
    warpx::subpixel::Mat3 m{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            m(i, j) = amrex::Real(0.5) * (a[i][j] + a[j][i]);
        }
    }
    return m;
}

}  // namespace

int main(int argc, char** argv)
{
    std::string in  = (argc > 1) ? argv[1] : "p1_cases.csv";
    std::string out = (argc > 2) ? argv[2] : "p1_cpp.csv";

    std::ifstream fin(in);
    if (!fin) { std::fprintf(stderr, "cannot open %s\n", in.c_str()); return 1; }
    std::ofstream fout(out);
    if (!fout) { std::fprintf(stderr, "cannot write %s\n", out.c_str()); return 1; }
    fout << "case,chi_00,chi_01,chi_02,chi_10,chi_11,chi_12,chi_20,chi_21,chi_22\n";
    fout.precision(17);

    std::string line;
    std::getline(fin, line);  // header
    int ncase = 0;
    while (std::getline(fin, line)) {
        if (line.empty()) { continue; }
        auto v = split(line);
        if (v.size() < 23) { continue; }
        int const case_id = std::stoi(v[0]);
        warpx::subpixel::Mat3 const e1 = sym_from_row(v, 2);
        warpx::subpixel::Mat3 const e2 = sym_from_row(v, 11);
        amrex::Real const n[3] = { std::stod(v[20]), std::stod(v[21]), std::stod(v[22]) };
        amrex::Real const f = std::stod(v[23]);

        amrex::Real chi[9];
        warpx::subpixel::chi_inv_anisotropic(e1, e2, f, n, chi);

        fout << case_id;
        for (int i = 0; i < 9; ++i) { fout << ',' << chi[i]; }
        fout << '\n';
        ++ncase;
    }
    std::printf("wrote %s (%d cases)\n", out.c_str(), ncase);
    return 0;
}
