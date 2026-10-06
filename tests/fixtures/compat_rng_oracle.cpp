// Target-native conformance probe for the deliberately retained C rand/srand
// dependency. This is compiled only by the opt-in RNG oracle test; it is not
// part of the Rust library or binary build.

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "lib/tantan/LambdaCalculator.hh"
#include "stats/matrices/blosum62.h"

static constexpr double PI = 3.1415926535897932384626433832795;

static double standard_normal() {
  double r1 = 0.0;
  while (r1 == 0.0) r1 = static_cast<double>(rand()) / RAND_MAX;
  double r2 = 0.0;
  while (r2 == 0.0) r2 = static_cast<double>(rand()) / RAND_MAX;
  double v1 = -2.0 * log(r1);
  if (v1 < 0.0) v1 = 0.0;
  return sqrt(v1) * cos(2.0 * PI * r2);
}

static void print_random_vector(unsigned seed) {
  srand(seed);
  std::printf("rand-%u", seed);
  for (int i = 0; i < 100; ++i) std::printf(" %d", rand());
  std::putchar('\n');
}

static void print_normal_vector() {
  srand(12345);
  std::printf("standard-normal");
  for (int i = 0; i < 16; ++i) std::printf(" %.17g", standard_normal());
  std::putchar('\n');
}

static void print_tantan_vector() {
  int matrix[20][20];
  const int* rows[20];
  for (int i = 0; i < 20; ++i) {
    rows[i] = matrix[i];
    for (int j = 0; j < 20; ++j)
      matrix[i][j] = Stats::blosum62.scores[i * 26 + j];
  }

  srand(1);
  cbrc::LambdaCalculator calculator;
  calculator.calculate(rows, 20);
  std::printf("tantan %.17g", calculator.lambda());
  for (int i = 0; i < 20; ++i)
    std::printf(" %.17g", calculator.letterProbs1()[i]);
  for (int i = 0; i < 20; ++i)
    std::printf(" %.17g", calculator.letterProbs2()[i]);
  std::putchar('\n');
}

static void print_evaluator_samples() {
  // Exact call order in AlignmentEvaluer::initParametersWithErrors: twelve
  // standard-normal draws per row. Emit the first and last field so every
  // intervening draw affects the comparison.
  double lambdas[20];
  double taus[20];
  const double scale = sqrt(20.0);
  srand(12345);
  for (int i = 0; i < 20; ++i) {
    lambdas[i] = 0.267 + 0.001 * standard_normal() * scale;
    (void)(0.041 + 0.001 * standard_normal() * scale);  // K
    (void)(43.0 + 0.01 * standard_normal() * scale);   // sigma
    (void)(1.8 + 0.01 * standard_normal() * scale);    // alpha_I
    (void)(1.7 + 0.01 * standard_normal() * scale);    // alpha_J
    (void)(2.1 + 0.01 * standard_normal() * scale);    // a_I
    (void)(1.9 + 0.01 * standard_normal() * scale);    // a_J
    (void)(5.0 + 0.01 * standard_normal() * scale);    // b_I
    (void)(4.0 + 0.01 * standard_normal() * scale);    // b_J
    (void)(7.0 + 0.01 * standard_normal() * scale);    // beta_I
    (void)(6.0 + 0.01 * standard_normal() * scale);    // beta_J
    taus[i] = 8.0 + 0.01 * standard_normal() * scale;
  }
  std::printf("evaluator");
  for (double value : lambdas) std::printf(" %.17g", value);
  for (double value : taus) std::printf(" %.17g", value);
  std::putchar('\n');
}

int main() {
  print_random_vector(1);
  print_random_vector(12345);
  print_normal_vector();
  print_tantan_vector();
  print_evaluator_samples();
  return 0;
}
