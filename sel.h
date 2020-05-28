float calculateLocalD(int i,mesh m){
    float D,a,b,c,d;

    element e = m.getElement(i);

    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);

    a = n2.getX()-n1.getX(); b = n2.getY()-n1.getY();
    c = n3.getX()-n1.getX(); d = n3.getY()-n1.getY();
    D = a*d - b*c;

    return D;
}

float calculateMagnitude(float v1, float v2){
    return sqrt(pow(v1,2)+pow(v2,2));
}

float calculateLocalArea(int i,mesh m){
    float A,s,a,b,c;
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);

    a = calculateMagnitude(n2.getX()-n1.getX(),n2.getY()-n1.getY());
    b = calculateMagnitude(n3.getX()-n2.getX(),n3.getY()-n2.getY());
    c = calculateMagnitude(n3.getX()-n1.getX(),n3.getY()-n1.getY());
    s = (a+b+c)/2;

    A = sqrt(s*(s-a)*(s-b)*(s-c));
    return A;
}

void calculateLocalA(int i,Matrix &A,mesh m){
    
    element e = m.getElement(i);
    
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);

    A.at(0).at(0) = n3.getY()-n1.getY(); A.at(0).at(1) = n1.getY()-n2.getY();
    A.at(1).at(0) = n1.getX()-n3.getX(); A.at(1).at(1) = n2.getX()-n1.getX();
}

void calculateB(Matrix &B){
    B.at(0).at(0) = -1; B.at(0).at(1) = 1; B.at(0).at(2) = 0;
    B.at(1).at(0) = -1; B.at(1).at(1) = 0; B.at(1).at(2) = 1;
}

Matrix createLocalK(int element,mesh &m){
    float D,Ae,k = m.getParameter(THERMAL_CONDUCTIVITY);
    Matrix K,A,B,Bt,At;

    D = calculateLocalD(element,m);

    if(D == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }

    Ae = calculateLocalArea(element,m);

    zeroes(A,2);
    zeroes(B,2,3);
    calculateLocalA(element,A,m);
    calculateB(B);
    transpose(A,At);
    transpose(B,Bt);

    productRealMatrix(k*Ae/(D*D),productMatrixMatrix(Bt,productMatrixMatrix(At,productMatrixMatrix(A,B,2,2,3),2,2,3),3,2,3),K);

    return K;
}

float calculateLocalJ(int i,mesh m){
    float J,a,b,c,d;
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);

    a = n2.getX()-n1.getX(); b = n3.getX()-n1.getX();
    c = n2.getY()-n1.getY(); d = n3.getY()-n1.getY();
    
    J = a*d - b*c;

    return J;
}

Vector createLocalb(int element,mesh &m){
    Vector b;

    float Q = m.getParameter(HEAT_SOURCE),J;
    J = calculateLocalJ(element,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }

    b.push_back(Q*J*(1/6)); 
    b.push_back(Q*J*(1/6)); 
    b.push_back(Q*J*(1/6));

    return b;
}

void crearSistemasLocales(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs){
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        localKs.push_back(createLocalK(i,m));
        localbs.push_back(createLocalb(i,m));
    }
}

void assemblyK(element e,Matrix localK,Matrix &K){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;
    int index3 = e.getNode3() - 1;

    K.at(index1).at(index1) += localK.at(0).at(0);
    K.at(index1).at(index2) += localK.at(0).at(1);
    K.at(index1).at(index3) += localK.at(0).at(2);
    K.at(index2).at(index1) += localK.at(1).at(0);
    K.at(index2).at(index2) += localK.at(1).at(1);
    K.at(index2).at(index3) += localK.at(1).at(2);
    K.at(index3).at(index1) += localK.at(2).at(0);
    K.at(index3).at(index2) += localK.at(2).at(1);
    K.at(index3).at(index3) += localK.at(2).at(2);
}

void assemblyb(element e,Vector localb,Vector &b){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;
    int index3 = e.getNode3() - 1;

    b.at(index1) += localb.at(0);
    b.at(index2) += localb.at(1);
    b.at(index3) += localb.at(2);
}
