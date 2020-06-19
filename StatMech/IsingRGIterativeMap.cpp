#include <iostream>
#include <string>

double f(double, double), g(double, double);

int main() {
    double X, Y;
    double newX, newY;
    std::cout << "Initial X: ";
    std::cin >> X;
    std::cout << "Initial Y: ";
    std::cin >> Y;
    std::string flag;
    do {
        std::cout << "Iterate (X, Y) => (X', Y')?";
        getline(std::cin, flag);
        if (flag.empty()) {
            newX = f(X, Y);
            newY = g(X, Y);
            X = newX;
            Y = newY;
            std::cout << X << "    " << Y << "\n";
        }
        else {
            break;
        }
    } while (true);
    return 0;
}

double f(double X, double Y){
    return X * (1 + Y) * (1 + Y) / ((X + Y) * (1 + X * Y));
}

double g(double X, double Y){
    return Y * (X + Y) / (1 + X * Y);
}

