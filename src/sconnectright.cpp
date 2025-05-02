
#include "../include/sconnectright.h"

void SconnectRight(SMatrix& G, SMatrix& S)
{
	

	cx_mat Identity(S.S11.n_rows, S.S11.n_cols, fill::eye);
	cx_mat D = G.S12 * inv(Identity - S.S11 * G.S22);
	cx_mat F = S.S21 * inv(Identity - G.S22 * S.S11);


	G.S11 = G.S11 + D * S.S11 * G.S21;
	G.S12 = D * S.S12;
	G.S21 = F * G.S21;
	G.S22 = S.S22 + F * G.S22 * S.S12;
}

//SMatrix sRight(const SMatrix& A, const SMatrix& B)
//{
//	SMatrix AB;
//	cx_mat II = eye(size(A.S11));
//	AB.S11 = A.S11 + A.S12 * inv(II - B.S11 * A.S22) * B.S11 * A.S21;
//	AB.S12 =  A.S12 * inv(II-B.S11 * A.S22) * B.S12  ;
//	AB.S21 = B.S21* inv(II- A.S22 * B.S11) *   A.S21;
//	AB.S22 = B.S22 + B.S21 * inv(II - A.S22 * B.S11) *A.S22*B.S12;
//	return AB;
//}



