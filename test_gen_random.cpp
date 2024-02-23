#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

// Fonction pour générer une particule dans un domaine cubique discrétisé
std::vector<double> genererParticuleCubique(double Lx, double Ly, double Lz, int Nx, int Ny, int Nz) {
    // Initialisation de rand()
    srand(static_cast<unsigned int>(time(nullptr)));

    // Génération aléatoire des coordonnées de la particule
    double x = static_cast<double>(rand()) / RAND_MAX * Lx;
    double y = static_cast<double>(rand()) / RAND_MAX * Ly;
    double z = static_cast<double>(rand()) / RAND_MAX * Lz;

    // Discrétisation des coordonnées
    int indexX = static_cast<int>((x / Lx) * Nx);
    int indexY = static_cast<int>((y / Ly) * Ny);
    int indexZ = static_cast<int>((z / Lz) * Nz);

    // Affichage des coordonnées discrètes de la particule
    std::cout << "Coordonnées discrètes de la particule : (" << indexX << ", " << indexY << ", " << indexZ << ")\n";

    // Retourner les coordonnées de la particule
    return {x, y, z};
}

int main() {
    double Lx = 10.0; // Longueur du côté en x
    double Ly = 10.0; // Longueur du côté en y
    double Lz = 10.0; // Longueur du côté en z
    int Nx = 5; // Nombre d'éléments le long de l'axe x
    int Ny = 5; // Nombre d'éléments le long de l'axe y
    int Nz = 5; // Nombre d'éléments le long de l'axe z

    // Génération d'une particule dans le domaine cubique discrétisé
    std::vector<double> particule = genererParticuleCubique(Lx, Ly, Lz, Nx, Ny, Nz);

    // Affichage des coordonnées réelles de la particule
    std::cout << "Coordonnées réelles de la particule : (" << particule[0] << ", " << particule[1] << ", " << particule[2] << ")\n";

    return 0;
}