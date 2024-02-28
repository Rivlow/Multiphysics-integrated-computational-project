#include <iostream>
#include <vector>

int main() {
    const int m = 3; // Nombre de colonnes (première dimension fixe)

    std::vector<std::vector<int>> tableau = {
        {1, 2, 3},
        {4, 5},
        {6, 7, 8, 9}
    };

    // Imprimer chaque ligne du tableau
    for (size_t i = 0; i < tableau.size(); ++i) {
        std::cout << "Ligne " << i << ": ";
        for (size_t j = 0; j < tableau[i].size(); ++j) {
            std::cout << tableau[i][j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}