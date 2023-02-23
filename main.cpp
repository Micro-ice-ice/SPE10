#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char **argv) {

    string filename = argv[0];   // Name of the file
    ifstream newfile (filename);
    newfile.open(filename);
    if (newfile.is_open()){ //checking whether the file is open
        cout << "File is open";
        string line;
        while (newfile >> line) {
            cout << line << endl;
        }
        newfile.close(); //close the file object.
    } else {
        cout << "File isn't open";
    }

    return 0;

}
