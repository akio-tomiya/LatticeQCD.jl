#include <Grid/Grid.h>
#include <Grid/qcd/smearing/HISQSmearing.h>

#include <complex>
#include <iomanip>
#include <iostream>

using namespace Grid;

namespace {

LatticeGaugeFieldD deterministic_links(GridCartesian *grid) {
  LatticeGaugeFieldD gauge(grid);
  gauge = Zero();
  for (int t = 0; t < 4; ++t) {
    for (int z = 0; z < 4; ++z) {
      for (int y = 0; y < 4; ++y) {
        for (int x = 0; x < 4; ++x) {
          const Coordinate site({x, y, z, t});
          LorentzColourMatrixD value;
          value = Zero();
          const int coordinate = (x + 1) + 3 * (y + 1) +
                                 5 * (z + 1) + 7 * (t + 1);
          for (int mu = 0; mu < Nd; ++mu) {
            for (int column = 0; column < Nc; ++column) {
              for (int row = 0; row < Nc; ++row) {
                const double re = 0.05 * 0.013 *
                    (2 * (row + 1) - (column + 1) + coordinate +
                     3 * (mu + 1));
                const double im = 0.05 * 0.017 *
                    ((row + 1) + 2 * (column + 1) - coordinate +
                     (mu + 1));
                value(mu)()(row, column) =
                    ComplexD(re + (row == column ? 1.0 : 0.0), im);
              }
            }
          }
          pokeSite(value, gauge, site);
        }
      }
    }
  }
  return gauge;
}

void print_fingerprints(const char *stage, const LatticeGaugeFieldD &gauge) {
  for (int mu = 0; mu < Nd; ++mu) {
    const auto link = PeekIndex<LorentzIndex>(gauge, mu);
    double sum_re = 0.0;
    double sum_im = 0.0;
    double weighted_re = 0.0;
    double weighted_im = 0.0;
    double norm2_value = 0.0;
    long index = 0;
    for (int t = 0; t < 4; ++t) {
      for (int z = 0; z < 4; ++z) {
        for (int y = 0; y < 4; ++y) {
          for (int x = 0; x < 4; ++x) {
            ColourMatrixD matrix;
            peekSite(matrix, link, Coordinate({x, y, z, t}));
            for (int column = 0; column < Nc; ++column) {
              for (int row = 0; row < Nc; ++row) {
                ++index;
                const ComplexD value =
                    TensorRemove(matrix()()(row, column));
                sum_re += value.real();
                sum_im += value.imag();
                weighted_re += index * value.real();
                weighted_im += index * value.imag();
                norm2_value +=
                    value.real() * value.real() + value.imag() * value.imag();
              }
            }
          }
        }
      }
    }
    std::cout << std::setprecision(17)
              << "FINGERPRINT implementation=Grid stage=" << stage
              << " mu=" << mu + 1
              << " sum_re=" << sum_re
              << " sum_im=" << sum_im
              << " weighted_re=" << weighted_re
              << " weighted_im=" << weighted_im
              << " norm2=" << norm2_value << '\n';
  }
}

}  // namespace

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
  const Coordinate lattice({4, 4, 4, 4});
  const Coordinate mpi({1, 1, 1, 1});
  auto *grid = SpaceTimeGrid::makeFourDimGrid(
      lattice, GridDefaultSimd(Nd, vComplexD::Nsimd()), mpi);

  LatticeGaugeFieldD thin = deterministic_links(grid);
  LatticeGaugeFieldD level1(grid), unused_naik(grid), reunitarized(grid);
  LatticeGaugeFieldD fat(grid), long_links(grid);

  Smear_HISQ<PeriodicGimplD> first(
      grid, 1.0 / 8.0, 0.0, 1.0 / 16.0, 1.0 / 64.0,
      1.0 / 384.0, 0.0);
  first.smear(level1, unused_naik, thin);
  first.projectU3(reunitarized, level1);

  constexpr double naik_epsilon = -0.083;
  Smear_HISQ<PeriodicGimplD> second(
      grid, 1.0 + naik_epsilon / 8.0, 1.0, 1.0 / 16.0,
      1.0 / 64.0, 1.0 / 384.0, -1.0 / 8.0);
  second.smear(fat, long_links, reunitarized);

  if (grid->IsBoss()) {
    print_fingerprints("fat", fat);
    print_fingerprints("long", long_links);
  }

  delete grid;
  Grid_finalize();
  return 0;
}
