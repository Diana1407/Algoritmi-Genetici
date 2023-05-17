#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace std;

//ifstream f("date.in"); ///fisier pt functia de gradul3
ifstream f("date2.in"); ///fisier pt functia de gradul 2
ofstream g("evol.out");

int n, c[4], precizie, nr_etape, len, ind, p;
long double maxx, p_recombinare, p_mutatie, minn, a, b;
vector <long long> cromozomi, new_gen;
vector <long double> vals, probs, intervale;
vector <int> viz;
vector <int> recombinare, aux;

/// functie de afisare a numerelor in baza 2
/// se parcurge nr de la cel mai semnificativ bit al sau si daca bitul e setat afisam 1 altfel 0
void transf_binar(long long nr) {
    long long i;
    for (i = 1 << (len - 1); i > 0; i = i / 2) {
        if((nr & i) != 0) {
            g << "1";
        }
        else {
            g << "0";
        }
    }
}

/// functie de rotunjire cu p zecimale
long double rotunjire (long double nr, int p) {
    int value = (int)(nr * (long double)p);
    return (long double)value / (long double)p;
}

/// functie ce caluleaza functia de gradul 2
/// functia pe care vrem sa o maximizam
long double functie (long double x) {
    return ((long double)c[0] * x * x + c[1] * x + c[2]);
}

long double functie3(long double x)
{
    return ((long double)c[0] * x * x * x  + (long double)c[1] * x * x + c[2] * x + c[3]);
}

int cautare_binara(long double x) {
    int mijl, st = 0, dr = n;
    while (st <= dr) {
        mijl = (st + dr);
        if (intervale[mijl] <= x) {
            st = mijl + 1;
        } else if (intervale[mijl] > x) {
            dr = mijl - 1;
        }
    }
    return dr;
}

void cit()
{
    //cout << "Introduceti dimensiunea populatiei: ";
    f >> n;
    //cout<<n<<endl;
    //cout << "Introduceti domeniul de definitie al functiei: ";
    f >> a >> b;
    //cout<<a<<" "<<b<<endl;
    //cout << "Introduceti coeficientii functiei: ";
    f >> c[0] >> c[1] >> c[2];
    //f>> c[0] >> c[1] >> c[2] >> c[3];
    //cout<<c[0]<<" "<<c[1]<<" "<<c[2]<<" "<<c[3];
    //cout << "Introduceti precizia cu care se lucreaza: ";
    f >> precizie;
    //cout << "Introduceti probabilitatea de recombinare: ";
    f >> p_recombinare;
    //cout << "Introduceti probabilitatea de mutatie: ";
    f >> p_mutatie;
    //cout << "Introduceti numarul de etape al algoritmului: ";
    f >> nr_etape;
    cout<<"Citire efectuata\n";
}

int main() {

    cit();

    /// formula din curs
    p = pow(10, precizie);
    len = ceil(log2((b - a) * p)); /// lungimea cromozomului
    long double d = ((long double)b - a) / (pow(2, len) - 1);


    /// la primul pas alegem indivizii populatiei random
    for (int i = 0; i < n; i++)
    {
        //int max_value = pow(2, len); echivalent cu 1 << len
        //cromozomi.push_back(rand() % max_value);
        cromozomi.push_back(rand() % (1 << len));
    }

    for (int gen = 0; gen < nr_etape; gen++) {
        if (gen == 0)
            g << "Populatia initiala:\n";

        long double total_fitness = 0;
        vals.clear();
        probs.clear();
        intervale.clear();
        viz.clear();
        recombinare.clear();

        long long ch_fit;
        for (int i = 0; i < n; i++) {
            long long y = cromozomi[i];
            long double val_y = d * y + a;  /// valoarea codificata din interval - translatie liniara
            vals.push_back(val_y);

            if (gen == 0) {
                g << i + 1 << ": ";
                transf_binar(y);
                g << " x= " << rotunjire(vals[i], p) << " f= " << functie(rotunjire(vals[i], p)) << '\n';
                //g << " x= " << rotunjire(vals[i], p) << " f= " << functie3(rotunjire(vals[i], p)) << '\n';
            }

            ///calculam fitness-ul total si updatam maximul si indicele maximului
            total_fitness += functie(vals[i]);
            //total_fitness += functie3(vals[i]);

            if (i == 0) {
                maxx = functie(rotunjire(vals[i], p));
                //maxx = functie3(rotunjire(vals[i], p));
                ind = i;
            }
            else if (functie(rotunjire(vals[i], p)) > maxx){
            //else if (functie3(rotunjire(vals[i], p)) > maxx) {
                //maxx = functie(rotunjire(vals[i], p));
                maxx = functie3(rotunjire(vals[i], p));
                ind = i;
            }
        }

        long double max_fit = maxx; /// cel mai mare scor fitness
        ch_fit = cromozomi[ind]; /// cel mai fit individ din generatie

        if (gen == 0) {
            g << "\nProbabilitati selectie:\n";
        }
        for (int i = 0; i < n; i++) {
            long double prob_i = functie(vals[i]) / total_fitness;
            //long double prob_i = functie3(vals[i]) / total_fitness;
            probs.push_back(prob_i);
            if (gen == 0) {
                g << "cromozom " << i + 1 << " probabilitate " << prob_i << '\n';
            }
        }

        if (gen == 0) {
            g << "\nIntervale probabilitati selectie:\n";
        }
        intervale.push_back(0);
        long double sum_interv = probs[0];
        intervale.push_back(sum_interv);

        if (gen == 0) {
            g << 0 << " " << sum_interv << " ";
        }
        for (int i = 1; i < n; i++) {
            sum_interv += probs[i];
            intervale.push_back(sum_interv);
            if (gen == 0) {
                g << sum_interv << ' ';
            }
        }

        /// selectia proportionala - metoda ruletei
        if (gen == 0)
            g<< '\n';
        for (int i = 0; i < n; i++) {
            /// generam un numar aleator din intervalul [0, 1) si ne asiguram ca este in acel interval
            /// prin impartirea la RAND_MAX care e cel mai mare nr random ce poate fi generat
            //long double u = static_cast <long double> (rand()) / static_cast <long double> (RAND_MAX);
            long double u = ((long double) rand())/ ((long double) RAND_MAX);
            int cromozom = cautare_binara(u);
            if (gen == 0) {
                g << "u= " << u << " selectam cromozomul " << cromozom + 1 << '\n';
            }
            viz.push_back(cromozom);
        }

        if (gen == 0) {
            g << "\nDupa selectie:\n";
        }
        for (int i = 0; i < n; i++) {
            if (gen == 0) {
                g << i + 1 << ": ";
                transf_binar(cromozomi[viz[i]]);
                g << " x= " << vals[viz[i]] << " f= " << functie(vals[viz[i]]) << '\n';
                //g << " x= " << vals[viz[i]] << " f= " << functie3(vals[viz[i]]) << '\n';
            }
            new_gen.push_back(cromozomi[viz[i]]);
        }
        cromozomi = new_gen;

        /// selectam cromozomii care vor fi recombinati
        if (gen == 0) {
            g << "\nProbabilitatea de incrucisare: " << p_recombinare << '\n';
        }
        for (int i = 0; i < n; i++) {
            //long double u = static_cast <long double> (rand()) / static_cast <long double> (RAND_MAX);
            long double u = ((long double) rand())/ ((long double) RAND_MAX);
            if (u < p_recombinare) {
                recombinare.push_back(i);
                if (gen == 0) {
                    g << i + 1 << ": ";
                    transf_binar(cromozomi[i]);
                    g << " u= " << u << " < " << p_recombinare << " participa\n";
                }
            }
            else {
                if (gen == 0) {
                    g << i + 1 << ": ";
                    transf_binar(cromozomi[i]);
                    g << " u= " << u << "\n";
                }
            }
        }

        /// recombinam cate 2 cromozomi random
        while (recombinare.size() > 1) {
            int x = rand() % recombinare.size();
            int y = rand() % recombinare.size();
            if (x == y) continue;
            if (gen == 0) {
                g << "Recombinare dintre cromozomul " << recombinare[x] + 1 << " cu cromozomul " << recombinare[y] + 1 << ":\n";
            }
            int taie = rand() % len;
            if (gen == 0) {
                transf_binar(cromozomi[recombinare[x]]);
                g << ' ';
                transf_binar(cromozomi[recombinare[y]]);
                g << " punct " << taie << "\nRezultat ";

            }
            long long masca = (1 << (len - taie + 1)) -1;
            long long masca_x = cromozomi[recombinare[x]] & masca;
            long long masca_y = cromozomi[recombinare[y]] & masca;
            ///taiem primii len - taie biti
            cromozomi[recombinare[x]] >>= (len - taie);
            ///facem loc pentru a putea dupa imbina
            cromozomi[recombinare[x]] <<= (len - taie);
            ///facem or pentru a adauga corect bitii din masca, adica pe pozitia corecta
            cromozomi[recombinare[x]] |= masca_y;
            cromozomi[recombinare[y]] >>= (len - taie);
            cromozomi[recombinare[y]] <<= (len - taie);
            cromozomi[recombinare[y]] |= masca_x;

            if (gen == 0) {
                transf_binar(cromozomi[recombinare[x]]);
                g << ' ';
                transf_binar(cromozomi[recombinare[y]]);
                g << '\n';
            }

            /// elimin cromozomii abia recombinati din vectorul de recombinare
            aux.clear();
            for (int i = 0; i < recombinare.size(); i++) {
                if (i != x && i != y) {
                    aux.push_back(recombinare[i]);
                }
            }
            recombinare = aux;
        }

        if (gen == 0) {
            g << "\nDupa recombinare:\n";
        }

        for (int i = 0; i < n; i++) {

            /// recalculez valoarea codificata din interval
            long long y = cromozomi[i];
            long double val_y = d * y + a;
            vals[i] = val_y;

            if (gen == 0) {
                g << i + 1 << ": ";
                transf_binar(y);
                g << " x= " << vals[i] << " f= " << rotunjire(functie(vals[i]), p) << '\n';
                //g << " x= " << vals[i] << " f= " << rotunjire(functie3(vals[i]), p) << '\n';
            }
        }

        if (gen == 0) {
            g << "\nProbabilitatea de mutatie " << p_mutatie << '\n';
            g << "Au fost modificati cromozomii:\n";
        }
        for (int i = 0; i < n; i++) {
            //long double u = static_cast <long double> (rand()) / static_cast <long double> (RAND_MAX);
            long double u = ((long double) rand())/ ((long double) RAND_MAX);
            if (u < p_mutatie) {
                int poz = rand() % len;
                long long masca = (1 << (len - poz + 1)) - 1;
                masca &= cromozomi[i]; /// aplicam masca pe cromozom pt a pastra intacta partea pe care nu vrem sa o schimbam
                cromozomi[i] >>= (len - poz); /// deplasam bitii la drepata
                if (cromozomi[i] % 2 == 0) /// verificam paritatea celui mai semnificativ bit pt a asigura schimbarea
                    cromozomi[i]++;
                else
                    cromozomi[i]--;
                cromozomi[i] <<= (len - poz); ///reintroducem bitii eliminatii
                cromozomi[i] |= masca; /// facem or cu masca pt a recombina partea schimbata cu cea neschimbata

                if (gen == 0) {
                    g << i + 1 << '\n';
                }
            }
        }

        if (gen == 0) {
            g << "\nDupa mutatie:\n";
        }

        minn = 0;
        maxx = 0;
        for (int i = 0; i < n; i++) {

            /// recalculez valoarea codificata din interval
            long long y = cromozomi[i];
            long double val_y = d * y + a;
            vals[i] = val_y;
            long double fct = functie(vals[i]);
            //long double fct = functie3(vals[i]);
            if (gen == 0) {
                g << i + 1 << ": ";
                transf_binar(y);
                g << " x= " << vals[i] << " f= " << rotunjire(fct, p) << '\n';
            }
            if (i == 0) {
                maxx = fct;
                minn = fct;
            }
            else {
                /// retin maximul de fitness
                if (fct > maxx) {
                    maxx = fct;
                }
                /// retin cel mai nefit cromozom
                if (fct < minn) {
                    minn = fct;
                    ind = i;
                }
            }
        }

        /// selectia elitista
        /// schimb cel mai nefit cromozom cu cel mai fit din prima generatie inainte de schimbari
        cromozomi[ind] = ch_fit;


        if (max_fit > maxx) {
            maxx = max_fit;
        }

        if (gen == 0) {
            g << "\nEvolutia maximului:\n";
        }

        g << maxx << '\n';

    }

    return 0;
}
