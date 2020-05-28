#include <fstream>
#include "string.h"

void obtenerDatos(istream &file,int nlines,int n,int mode,item* item_list){
    string line;
    file >> line;
    if(nlines==DOUBLELINE) file >> line;

    for(int i=0;i<n;i++){
        switch(mode){
        case INT_FLOAT:
            int e0; float condition_value;
            file >> e0 >> condition_value;
            item_list[i].setValues(0,0,0,0,0,e0,condition_value,0);
            break;
        case INT_FLOAT_FLOAT_FLOAT:
            int e; float node1, node2, node3;
            file >> e >> node1 >> node2 >> node3;
            item_list[i].setValues(e,node1,node2,node3,0,0,0,0);
            break;
        case INT_INT_INT_INT_INT:
            int element,e1,e2,e3,e4;
            file >> element >> e1 >> e2 >> e3 >> e4;
            item_list[i].setValues(element,0,0,0,e1,e2,e3,e4);
            break;
        }
    }
}

void correctConditions(int n,condition *list,int *indices){
    for(int i=0;i<n;i++)
        indices[i] = list[i].getNode1();

    for(int i=0;i<n-1;i++){
        int pivot = list[i].getNode1();
        for(int j=i;j<n;j++)
            if(list[j].getNode1()>pivot)
                list[j].setNode1(list[j].getNode1()-1);
    }
}
    
void addExtension(char *newfilename,char *filename, string extension){
    int ori_length = strlen(filename);
    int ext_length = extension.length();
    int i;
    for(i=0;i<ori_length;i++)
        newfilename[i] = filename[i];
    for(i=0;i<ext_length;i++)
        newfilename[ori_length+i] = extension[i];
    newfilename[ori_length+i] = '\0';
}

void leerMallayCondiciones(mesh &m,char *filename){
    char inputfilename[150];
    ifstream file;
    float k,Q_x, Q_y, Q_z;
    int nnodes,neltos,ndirich,nneu;

    addExtension(inputfilename,filename,".dat");
    file.open(inputfilename);

    file >> k >> Q_x >> Q_y >> Q_z;

    file >> nnodes >> neltos >> ndirich >> nneu;

    m.setParameters(k,Q_x, Q_y, Q_z);
    m.setSizes(nnodes,neltos,ndirich,nneu);
    m.createData();

    obtenerDatos(file,SINGLELINE,nnodes,INT_FLOAT_FLOAT_FLOAT,m.getNodes());
    obtenerDatos(file,DOUBLELINE,neltos,INT_INT_INT_INT_INT,m.getElements());
    obtenerDatos(file,DOUBLELINE,ndirich,INT_FLOAT,m.getDirichlet());
    obtenerDatos(file,DOUBLELINE,nneu,INT_FLOAT,m.getNeumann());

    file.close();

    correctConditions(ndirich,m.getDirichlet(),m.getDirichletIndices());
}

bool findIndex(int v, int s, int *arr){
    for(int i=0;i<s;i++)
        if(arr[i]==v) return true;
    return false;
}

void writeResults(mesh m,Vector T,char *filename){
    char outputfilename[150];
    int *dirich_indices = m.getDirichletIndices();
    condition *dirich = m.getDirichlet();
    ofstream file;

    addExtension(outputfilename,filename,".post.res");
    file.open(outputfilename);

    file << "GiD Post Results File 1.0\n";
    file << "Result \"Temperature\" \"Load Case 1\" 1 Scalar OnNodes\nComponentNames \"T\"\nValues\n";

    int Tpos = 0;
    int Dpos = 0;
    int n = m.getSize(NODES);
    for(int i=0;i<n;i++){
        if(findIndex(i+1,n,dirich_indices)){
            file << i+1 << " " << dirich[Dpos].getValue() << "\n";
            Dpos++;
        }else{
            file << i+1 << " " << T.at(Tpos) << "\n";
            Tpos++;
        }
    }

    file << "End values\n";

    file.close();
}
