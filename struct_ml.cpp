#include "struct_ml.h"
#include "struct_filion.h"
#include "struct_base.h"
#include "dyn_bb.h"
#include "defs.h"
#include "pbc.h"
#include <set>

// The comparison function for sorting the set by increasing order of its pair's
// second value. If the second value is equal, order by the pair's first value
struct comp
{
    template<typename T>
    bool operator()(const T &l, const T &r) const
    {
        if (l.second != r.second) {
            return l.second < r.second;
        }
 
        return l.first < r.first;
    }
};

void eval_struct_ml(){

    std::cout << "EVAL ML 1 " << std::endl; 

    int neigh_max_loc = 250;
    int ** neighbors = imatrix(0,NS*N-1,0,neigh_max_loc-1);
    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) {
            for (int d=0; d<neigh_max_loc;d++) {
                neighbors[i+s*N][d] = -1;
            }
        }
    }

    // calc neighbors
    double rcut = 6.0;
    double rcut2 = rcut*rcut;
    findneighbors(rcut2, neighbors);

    std::cout << "EVAL ML 2 " << std::endl; 

    // calc epot
    double dr, dx;
    double * epot = dvector(0,NS*N-1);
    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) {
            //std::cout << i << std::endl;
            int n0 = 0;
            while (neighbors[i+s*N][n0] != -1 ) {
                if (n0 > neigh_max_loc -2) std::cout << n0 << std::endl;
                int j= neighbors[i+s*N][n0];
                dr = 0;
                for (int d=0; d<dim;d++) {
                    dx = xyz_data[i+s*N][d] - xyz_data[j+s*N][d];
                    apply_pbc(dx);
                    dr += dx*dx;
                }
                // calc epot
                if (dr < RC2LJ) epot[i+s*N] += 0.5*calc_epot(i+s*N, j+s*N, dr);

                n0 ++;
            }
        }
    }

    std::cout << "EVAL ML 3 " << std::endl; 

    // create classifier ml
    double ** ml_classifiers = dmatrix(0,NS*N-1,0,3*neigh_max_loc - 1);
    for (int s=0; s<NS;s++) { // loop over structures
        for (int i=0; i<N;i++) {
            int n0 = 0;

            // create map
            std::map<int, double> map;

            while (neighbors[i+s*N][n0] != -1 ) {
                int j= neighbors[i+s*N][n0];
                dr = 0;
                for (int d=0; d<dim;d++) {
                    dx = xyz_data[i+s*N][d] - xyz_data[j+s*N][d];
                    apply_pbc(dx);
                    dr += dx*dx;
                }
                // include distance to map
                map[j] = dr;

                n0 ++;
            }

            // sort map by distance
            std::set<std::pair<int, double>, comp> set(map.begin(), map.end());

            // save result
            int counter = 0;
            for (auto const &pair: set) {
                int j = pair.first;
                double dist = pair.second;
                ml_classifiers[i+s*N][3*counter] = type_data[j+s*N];
                ml_classifiers[i+s*N][3*counter+1] = dist;
                ml_classifiers[i+s*N][3*counter+2] = epot[j+s*N];
                counter++;
            }

        }
    }

    std::cout << "EVAL ML 4 " << std::endl; 

    // print classifiers
    int print_max = 10;
    for (int type=0; type<NTYPE; type++) {
        QString pathPred = QString::fromStdString(folderOut);
        pathPred.append("struct_ml_type");
        pathPred.append(QString("%1.csv").arg(type));
        QFile outfilePred(pathPred);   // input file with xyz
        outfilePred.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream outPred(&outfilePred);
        //write header
        for (int typeJ=0; typeJ<NTYPE; typeJ++) {
            for (int print=0; print< print_max; print++) {
                outPred << "NT" << typeJ << ":" << print << ":DIST," << "NT" << typeJ << ":" << print << ":EPOT,";
            }
        }
        outPred << "\n";

        // write body
        
        // particle data
        for (int i=0; i<N*NS; i++) {
            if (type_data[i] == type) {
                for (int typeJ=0; typeJ<NTYPE; typeJ++) {
                    int count =0;
                    for (int print=0; print< neigh_max_loc; print++) {
                        if (count < print_max && typeJ == ml_classifiers[i][3*print]) {
                            outPred << ml_classifiers[i][3*print+1] << "," << ml_classifiers[i][3*print+2] << ",";
                            //std::cout << print << " " << count << std::endl;
                            count++;
                        }
                    }
                }
                outPred << "\n";
            }
        }
        outfilePred.close();
    }


}