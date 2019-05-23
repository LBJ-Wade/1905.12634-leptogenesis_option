#include "mr.hpp"

using namespace mr;

void output(SMCouplings av)
{
  std::cout << " g_1 = " << sqrt(av[couplings::g1]) * 4 * Pi
            << " g_2 = " << sqrt(av[couplings::g2]) * 4 * Pi
            << " g_3 = " << sqrt(av[couplings::gs]) * 4 * Pi
            << " y_t = " << sqrt(av[couplings::yt]) * 4 * Pi
            << " y_b = " << sqrt(av[couplings::yb]) * 4 * Pi
            << " y_tau = " << sqrt(av[couplings::ytau]) * 4 * Pi
            << " lam = " << av[couplings::lam] * 4 * Pi * 4 * Pi
            << " m = " << av[couplings::mphi] << std::endl;
}

int main(int argc, char *argv[])
{
  loglevel = logINFO;
  std::cout << std::setprecision(3);

  double
      alphas = 0.1185,
      gf = 1.1663787e-5;
  int as_loop_order = 2;
  int orders[3] = {0, 3, 31};

  double qs[2] = {173.2, 1e7};

  OSinput sminput(4.18, 80.387, 91.1875, 125.09, 173.2);
  AlphaS as(sminput, alphas, as_loop_order);

  for (int order = 0; order <= 2; order++)
  {
    P2MS<AlphaSolve> param_at_mt(sminput, gf, as(sminput.Mt()), sminput.Mt(), orders[order]);
    printf("%5d\t%.5Lf\t%.4Lf\t%.4Lf\t%.4Lf\t%.7Lf\t%.4Lf\t%.4Lf\n",
           order,
           param_at_mt.alam() * 4 * Pi * 4 * Pi,
           param_at_mt.mphi(),
           sqrt(param_at_mt.a1()) * 4 * Pi,
           sqrt(param_at_mt.a2()) * 4 * Pi,
           sqrt(param_at_mt.as()) * 4 * Pi,
           sqrt(param_at_mt.at()) * 4 * Pi,
           sqrt(param_at_mt.ab()) * 4 * Pi);
    if (order == 1)
    {
      ParametersSM<2, 2, 2, 2, 0, -1, 2, 2, 0> avP2MS(param_at_mt);
      for (size_t qi = 0; qi <= 1; qi++)
      {
        SMCouplings av = avP2MS(qs[qi] * qs[qi]);
        std::cout << "Q=" << qs[qi];
        output(av);
      }
    }
    if (order == 2)
    {
      ParametersSM<3, 3, 3, 3, 0, -1, 3, 3, 0> avP2MS(param_at_mt);
      for (size_t qi = 0; qi <= 1; qi++)
      {
        SMCouplings av = avP2MS(qs[qi] * qs[qi]);
        std::cout << "Q=" << qs[qi];
        output(av);
      }
    }
  }
}
