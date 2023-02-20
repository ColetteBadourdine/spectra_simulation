#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void distance_c(Rcpp::NumericMatrix reflectances, Rcpp::NumericMatrix distances) {
  for (int i=0; i < (reflectances.nrow()-1); i++) {
    for (int j=i+1; j < reflectances.nrow(); j++) {
      // Calcul de la distance
      distances(i,j) = sqrt(sum(pow(reflectances(i, _)-reflectances(j, _), 2)));
      // Symétrisation
      distances(j,i) = distances(i,j);
    }
  }
}


/*** R
# Appel de la fonction
#just a reminder to how to call the funtion

points_df <- data.frame(id = 1:4, x = rep(1:2, 2),
                        y = rep(1:2, each = 2), R = c(0, 255, 0, 0), V = c(0,
                                                                           0, 255, 0), B = c(0, 0, 0, 255))
points_n <- nrow(points_df)
points_d <- matrix(0, nrow = points_n, ncol = points_n)
distance_c(
  # Le premier argument est une matrice qui contient une ligne par pixel et tous les canaux en colonne
  as.matrix(points_df[, c("R", "V", "B")]),
  # Le deuxième argument est la matrice de résultat, qui doit être crée avant et avoir la bonne taille
  # Si la taille n'est pas bonne, R crashe.
  points_d
)
points_d
*/
