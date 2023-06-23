#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

void raytracer(const char *config, int show);

int main(int argc, const char *argv[]) {
  int show = 0;
  char *config = NULL;

  for (int i = 0; i < argc; ++i) {
    char *arg = (char *)argv[i];

    if ('-' == *arg) {
      ++arg;

      if ('h' == *arg) {
        printf(
          "usage: %s [-s] [config-file]\n"
          "options:\n"
          "  -h  display this page\n"
          "  -s  display traces\n",
          argv[0]
        );

        return EXIT_SUCCESS;
      } else if ('s' == *arg) {
        show = 1;
      } else {
        printf("error: '-%s' is not a valid option.\n", arg);
      }
    } else {
      config = arg;
    }
  }

  raytracer(config, show);

  return EXIT_SUCCESS;
}

void raytracer(const char *config, int show) {
  srand(time(NULL));

  const int MAXSC = 1000;

  long
    nrays = 100000 ,
    nscin = 3      ,
    nacts = 0      ,
    nhits = 0      ,
    n     = 10     ,
    nshit[MAXSC]   ;

  if (show) {
    nrays = 50;
  }

  if (nscin > MAXSC) {
    nscin = MAXSC;
  }

  for (long i = 0; i < nscin; ++i) {
    nshit[i] = 0;
  }

  const double PI2 = M_PI / 2;

  double
    l          = 0.1                 ,
    g          = 0.003               ,
    w          = 0.03                ,
    h          = 0.075               ,
    x_max      = n * l + (n - 1) * l ,
    x_mtor     = x_max / RAND_MAX    ,
    x                                ,
    theta                            ,
    theta_mtor = M_PI / RAND_MAX     ,
    theta_cos                        ,
    theta_prob                       ,
    theta_mean = 0.0                 , // mean of ALL theta generates
    ent_p0                           ,
    esc_p0                           ,
    q_slot     = l + g               , // our unit for measurements
    q_ent                            , // distance of ent_p0 from 0 in q_slots
    q_esc                            , // distance of esc_p0 from 0 in q_slots
    q_ps_ent                         , // position of ent_p0 in the q_slot
    q_ps_esc                         , // position of esc_p0 in the q_slot
    trace                            ,
    trace_aux                        , // used if entrance and escape are in different cells
    tr_act     = w / 10              ,
    act_eps    = 0.05                ;

  if (config) {
    FILE *fp = fopen(config, "r");

    if (!fp) {
      printf(
        "error: cannot open '%s'.\n"
        "note: RayTracer uses default data.\n",
        config
      );
    } else {
      fscanf(
        fp, "%ld%ld%ld%g%g%g%g%g",
        &nrays, &nscin, &n, &l, &g, &w, &h, &act_eps
      );
      fclose(fp);
    }
  }

  for(long i = 0; i < nrays; ++i) {

    // randomize x (linear)

    x = (double)rand() * x_mtor;

    // randomize theta (cos^2(x) distribution)

    do {
      theta_prob = (double)rand() / RAND_MAX;
      theta = (double)rand() * theta_mtor - PI2;
      theta_cos = cos(theta);
    } while (theta_prob > theta_cos * theta_cos);

    theta_mean += theta;

    // tracing...
    if (show) {
      printf("ray %ld - (%.2g, %.2g PI)\n", i + 1, x, theta / M_PI);
    }

    nhits = 0;
    for(int j = nscin; j > 0; --j) {
      // compute the points in the x-O reference system

      ent_p0 = x + ((j + 1) * h + (j + 1) * w) * tan(theta);
      esc_p0 = x + ((j + 1) * h + j * w) * tan(theta);

      // from x-O to 0-O reference system
      // ent_p0 = x + ent_px;
      // esc_p0 = x + esc_px;

      trace_aux = 0.0;

      // check if the entrance is inside or outside the scint.

      if (ent_p0 < 0 || ent_p0 > x_max) {
        // OR ent_px < 0 AND |ent_px| > x
        // OR ent_px > 0 AND ent_px > x_max - x

        // check oblique entrance

        if (esc_p0 < 0 || esc_p0 > x_max) {
          // no trace
          trace = 0.0;
        } else {

          if (ent_p0 < esc_p0) {
            trace = esc_p0;
          } else {
            trace = x_max - esc_p0;
          }

          trace /= sin(theta);
        }
      } else {
        // OR (ent_px <  0 AND |ent_px| < x)
        // OR ent_px >= 0

        q_ent = ent_p0 / q_slot;
        q_esc = esc_p0 / q_slot;

        q_ps_ent = (q_ent - (long)q_ent) * q_slot;
        q_ps_esc = (q_esc - (long)q_esc) * q_slot;

        if (q_ps_ent > l) {
          // entrance is in the gap (SURE)

          if (q_ps_esc > l) {
            if ((long)q_ent != (long)q_esc) {
              trace = l / sin(theta);
            } else {
              // no trace
              trace = 0.0;
            }
          } else {
            // escape in the cell (SURE)

            if (ent_p0 < esc_p0) {
              trace = q_ps_esc;
            } else {
              trace = l - q_ps_esc;
            }

            trace /= sin(theta);
          }
        } else {
          // entrace is in the cell (SURE)

          if (q_ps_esc > l) {
            // escape is in the gap (SURE)

            if (ent_p0 < esc_p0) {
              trace = l - q_ps_ent;
            } else {
              trace = q_ps_ent;
            }

            trace /= sin(theta);
          } else {
            // escape is in the cell (SURE)

            if ((long)q_ent != (long)q_esc) {
              // different cells

              if (ent_p0 < esc_p0) {
                trace = q_ps_ent;
                trace_aux = l - q_ps_esc;
              } else {
                trace = l - q_ps_ent;
                trace_aux = q_ps_esc;
              }

              trace /= sin(theta);
              trace_aux /= sin(theta);

              if (trace_aux < 0) {
                trace_aux = -trace_aux;
              }
            } else {
              // same cell

              trace = ent_p0 - esc_p0;
              trace = sqrt(w * w + trace * trace);
            }
          }
        }
      }

      // check activation

      if (trace < 0) {
        trace = -trace;
      }

      // show the ray
      if (show) {
        printf("sci. %2ld | ", j);
        for (long k = 0; k < n; ++k) {
          double lb, ub;

          if (k) {
            lb = k * l + (k - 1) * g;
            ub = lb + l;
          } else {
            lb = 0;
            ub = l;
          }

#ifdef ANSI_CONSOLE
          if (ent_p0 >= lb && ent_p0 < ub) {
            if (trace >= tr_act || trace >= tr_act) {
              double e = (double)rand() / RAND_MAX;

              if (e >= act_eps) {
                ++nhits;
                ++nshit[j - 1];
                printf("\033[92m1\033[0m");
              } else {
                printf("\033[93m2\033[0m");
              }
            } else {
              printf("\033[91m3\033[0m");
            }
          } else {
            printf("\033[90m0\033[0m");
          }
#else
          if (ent_p0 >= lb && ent_p0 < ub) {
            if (trace >= tr_act || trace >= tr_act) {
              double e = (double)rand() / RAND_MAX;

              if (e >= act_eps) {
                ++nhits;
                ++nshit[j - 1];
                putchar('1');
              } else {
                putchar('2');
              }
            } else {
              putchar('3');
            }
          } else {
            putchar('0');
          }
#endif
        }
        putchar('\n');

      } else {
        if (trace > tr_act || trace > tr_act) {
          double e = (double)rand() / RAND_MAX;

          if (e >= act_eps) {
            ++nhits;
            ++nshit[j - 1];
          }
        }
      }
    }

    if (nhits >= 2) {
      ++nacts;
    }
  }

  printf(
    ">>> RESULTS\n"
    "... act. = %ld\n"
    "... ang. = %.2e\n"
    "... eff. = %.2g%%\n",
    nacts, theta_mean / nrays, (double)nacts / nrays * 100
  );

  for (long i = 0; i < nscin; ++i) {
    printf("... sci. %2ld | %10ld | eff. = %.2g%%\n", i + 1, nshit[i], (double)nshit[i] / nrays * 100);
  }
}
