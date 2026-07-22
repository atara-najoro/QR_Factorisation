def produitScalaire(A, B):
        m = len(A)
        s = 0
        for i in range(m):
            s += A[i] * B[i]
        return s

def sommLigne(A, B):
        m = len(A)
        somme = [0] * m
        for i in range(m):
            somme[i] = A[i] + B[i]
        return somme

def normeVecteur(A):
        m = len(A)
        s = 0
        for i in range(m):
            s += A[i] * A[i]
        return pow(s, 1 / 2)

def sousLign(A, B):
        m = len(A)
        Res = [0] * m
        for i in range(m):
            Res[i] = A[i] - B[i]
        return Res

def MultiplicationParScalaire(coef, orig):
        m = len(orig)
        prod = [0] * m
        for i in range(m):
            prod[i] = coef * orig[i]
        return prod

def divisionParScalaire(div, A):
        m = len(A)
        quotient = [0] * m
        for i in range(m):
            quotient[i] = A[i] / div
        return quotient

def TransposeMatrice(Original):
        m = len(Original)
        n = len(Original[0])
        tO = [[0] * m for _ in range(n)]
        for i in range(n):
            for j in range(m):
                tO[i][j] = Original[j][i]
        return tO

def produitMat(prim, sec):
        a = len(prim)
        b = len(prim[0])
        d = len(sec[0])
        pro = [[0] * d for _ in range(a)]
        for i in range(a):
            for k in range(d):
                s = 0 
                for j in range(b):
                    s += prim[i][j] * sec[j][k]
                pro[i][k] = s
        return pro
def produit_scalaire(u, v):
        return sum(x * y for x, y in zip(u, v))

def produitMatriceVecteur(M, v):
        return [produit_scalaire(ligne, v) for ligne in M]

def decompositionQR(original):
        m, n = len(original), len(original[0])
        V = [[0] * m for _ in range(n)]
        A = TransposeMatrice(original)
        tR = [[0] * n for _ in range(n)]
        tQ = [[0] * m for _ in range(n)]
        for i in range(n):
            if i == 0:
                nor = normeVecteur(A[i])
                tR[0][0] = nor
                V[i] = A[i]
                if nor != 0:
                    tQ[i] = divisionParScalaire(nor, V[i])
            else:
                Sn = [0] * m 
                # Orthogonalisation
                for j in range(i):
                    r_val = produitScalaire(tQ[j], A[i])
                    tR[j][i] = r_val
                    Neg = MultiplicationParScalaire(r_val, tQ[j])
                    Sn = sommLigne(Sn, Neg)
                V[i] = sousLign(A[i], Sn)
                nor = normeVecteur(V[i])
                tR[i][i] = nor
                if nor > 1e-10:
                    tQ[i] = divisionParScalaire(nor, V[i])
                else:
                    tQ[i] = [0] * m
        tPReduit= [row for row in tQ if any(abs(x) > 1e-10 for x in row)]
        Q = TransposeMatrice(tPReduit)
        R = tR
        return Q, R