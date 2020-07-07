#include <iostream>
#include <string>
#include <cmath>

double f(double X, double Y){
    return X * (1 + Y) * (1 + Y) / ((X + Y) * (1 + X * Y));
}

double g(double X, double Y){
    return Y * (X + Y) / (1 + X * Y);
}

int main(int argc, char *argv[]) {
    double X = 0., Y = 0.;
    double newX = 1., newY = 1.;

    for (int i=1; i<argc; i+=2 ){
        std::string flag = std::string(argv[i]);
        std::string value = std::string(argv[i+1]);
        if (flag == "-x"){
            newX = stod(value);
        } else if (flag == "-y"){
            newY = stod(value);
        }
    }

    std::cout << newX << "    " << newY << "\n";
    while ((std::abs(newX-X) > pow(10.,-8.)) || (std::abs(newY-Y) > pow(10.,-8.))) {
        X = newX;
        Y = newY;
        newX = f(X, Y);
        newY = g(X, Y);
        std::cout << newX << "    " << newY << "\n";
    }
    return 0;
}