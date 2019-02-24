#include "mr.hpp"

using namespace mr;

int main(int argc, char *argv[])
{
  loglevel = logINFO;

  /*
  OSinput table2(4.2, 80.387, 91.1876, 125.15, 173.34);
  double
      alphas_table2 = 0.1184,
      gf_table2 = 1.166378139e-5;
  */
  OSinput table2(4.18, 80.387, 91.1875, 125.09, 173.2);
  double
      alphas_table2 = 0.1185,
      gf_table2 = 1.1663787e-5;
  int as_loop_order = 2;

  int orders[3] = {0, 3, 31};

  printf("order\tlam\tm\t\tg1\tg2\tg3      \tyt\tyb\n");
  for (int order = 0; order <= 2; order++)
  {
    AlphaS as(table2, alphas_table2, as_loop_order);
    P2MS<AlphaGF> param_at_mt(table2, gf_table2, as(table2.Mt()), table2.Mt(), orders[order]);
    printf("%5d\t%.5Lf\t%.4Lf\t%.4Lf\t%.4Lf\t%.7Lf\t%.4Lf\t%.4Lf\n",
           order,
           param_at_mt.alam() * 4 * Pi * 4 * Pi,
           param_at_mt.mphi(),
           sqrt(param_at_mt.a1()) * 4 * Pi,
           sqrt(param_at_mt.a2()) * 4 * Pi,
           sqrt(param_at_mt.as()) * 4 * Pi,
           sqrt(param_at_mt.at()) * 4 * Pi,
           sqrt(param_at_mt.ab()) * 4 * Pi);
  }
}
