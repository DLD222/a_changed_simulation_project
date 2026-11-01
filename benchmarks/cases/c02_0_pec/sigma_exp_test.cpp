/* Compare sigma_exponential_matrices() with the numpy reference.
 *
 * Prints, for a few fixed cases, the 18 entries of (M1, M2):
 *   case, m1_00..m1_22, m2_00..m2_22
 * sigma_exp_reference.py reads the same cases and prints the difference.
 */
#include <cstdio>

#include "SubpixelTensor.H"

using warpx::subpixel::Mat3;

namespace {

Mat3 from_rows(double a, double b, double c,
               double d, double e, double f,
               double g, double h, double i)
{
    Mat3 m{};
    m(0,0) = a; m(0,1) = b; m(0,2) = c;
    m(1,0) = d; m(1,1) = e; m(1,2) = f;
    m(2,0) = g; m(2,1) = h; m(2,2) = i;
    return m;
}

void emit(int id, const Mat3& eps, const Mat3& sigma, double dt)
{
    Mat3 M1{}, M2{};
    warpx::subpixel::sigma_exponential_matrices(eps, sigma, dt, M1, M2);
    std::printf("%d", id);
    for (int i = 0; i < 9; ++i) { std::printf(",%.17g", M1.m[i]); }
    for (int i = 0; i < 9; ++i) { std::printf(",%.17g", M2.m[i]); }
    std::printf("\n");
}

}  // namespace

int main()
{
    double const dt = 1.0e-12;
    Mat3 const I  = from_rows(1,0,0, 0,1,0, 0,0,1);
    Mat3 const eps4 = from_rows(4,0,0, 0,4,0, 0,0,4);

    // 0: sigma = 0 -> M1 = I, M2 = dt eps^-1
    emit(0, eps4, from_rows(0,0,0, 0,0,0, 0,0,0), dt);

    // 1: scalar sigma
    emit(1, eps4, from_rows(2,0,0, 0,2,0, 0,0,2), dt);

    // 2: diagonal anisotropic sigma
    emit(2, eps4, from_rows(1,0,0, 0,4,0, 0,0,9), dt);

    // 3: rotated anisotropic sigma (eps still isotropic)
    emit(3, eps4, from_rows(5.8,1.4,0, 1.4,4.2,0, 0,0,2.0), dt);

    // 4: anisotropic eps and sigma (both non-diagonal)
    emit(4, from_rows(5,1,0.5, 1,4,0.2, 0.5,0.2,3),
            from_rows(3,0.7,0, 0.7,6,0.1, 0,0.1,2.5), dt);

    // 5: large sigma*dt (stiff limit)
    emit(5, eps4, from_rows(1e9,0,0, 0,1e9,0, 0,0,1e9), dt);
    return 0;
}
