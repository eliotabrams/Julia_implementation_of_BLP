for n = 1:N 
	denom[n] =     sum{
        exp(beta[1]
            -(alpha+piInc[1]*inc[n]+piAge[1]*age[n]+sigma[1]*v[n,1])*p[h]
            +sum{(beta[k]+piInc[k]*inc[n]+piAge[k]*age[n]+sigma[k]*v[n,k])*x[h,k],k=2:K}
            +xi[h]
            )
    , h=1:J}

	for j = 1:J 																				# defining the tau as in the MPEC appendix
	tau[n,j] = exp(beta[1]
                -(alpha+piInc[1]*inc[n]+piAge[1]*age[n]+sigma[1]*v[n,1])*p[j]
                +sum{(beta[k]+piInc[k]*inc[n]+piAge[k]*age[n]+sigma[k]*v[n,k])*x[j,k],k=2:K}
                +xi[j]
            )/denom[n] 
	end
end

# theta = (alpha, beta, piInc, piAge, sigma). Will take derivative w.r.t. (theta, xi) => 1 + K + K + K + K + J variables

num_var = 1 + 4*K + J + L #number of variables oj jacobian: (theta, xi, g). Jacobian form: [ds/dtheta  ds/dxi  0 \\ 0 -Z' eye(g) ]

################################
# "Easy to understand version" #
################################

for row = 1:J # each row is the derivative of s_j wrt (theta,xi);

	values[1 + (row-1)*num_var] = (1/N)*sum{ tau[n,row]*(p[row] - sum{ tau[n,j]*p[j] , j=1:J } ) , n=1:N} #derivatives of s_j w.r.t. alpha
	values[2 + (row-1)*num_var] = 0 # derivative w.r.t. beta_0

	for column_beta = 2:K 
		values[1 + column_beta + (row-1)*num_var] = (1/N)*sum{ tau[n, row]*( x[row, column_beta] - sum{ tau[n,j]*x[j,column_beta] , j=1:J} , n=1:N) } #w.r.t to beta_2, beta_3, beta_4
	end

	values[2 + K + (row-1)*num_var] = (1/N)*sum{ tau[n,row]*( p[row] - sum{ tau[n,j]*p[j] , j=1:J } )*inc[n] , n=1:N } #piInc(1) => section that is multiplied by p

	for column_piInc = 2:K 
		values[1 + K + column_piInc + (row-1)*num_var] = (1/N)*sum{ tau[n, row]*(x[row,column_piInc] - sum{ tau[n,j]*x[j,column_piInc])*inc[n] , j=1:J } , n=1:N ) }
		#piInc(2), piInc(3) and piInc(4), the section multiplied by x
	end

	values[2 + 2*K + (row-1)*num_var] = (1/N)*sum{ tau[n,row]*( p[row] - sum{ tau[n,j]*p[j] , j=1:J } )*age[n] , n=1:N } #piAge(1) => section that is multiplied by p

	for column_piAge = 2:K 
		values[1 + 2*K + column_piAge + (row-1)*num_var] = (1/N)*sum{ tau[n, row]*(x[row,column_piAge] - sum{ tau[n,j]*x[j,column_piAge])*age[n] , j=1:J } , n=1:N ) }
		#piAge(2), piAge(3) and piAge(4), the section multiplied by x
	end

	values[2 + 3*K + (row-1)*num_var] = (1/N)*sum{ tau[n,row]*( p[row] - sum{ tau[n,j]*p[j] , j=1:J } )*v[n,1] , n=1:N } # w.r.t sigma_alpha

	for column_sigma_beta = 2:K 
		values[1 + 3*K + column_sigma_beta + (row-1)*num_var] = (1/N)*sum{ tau[n, row]*(x[row,column_sigma_beta] - sum{ tau[n,j]*x[j,column_sigma_beta])*v[n,column_sigma_beta] , j=1:J } , n=1:N ) } 
		#sigma_beta_2, sigma_beta_3, sigma_beta_4
	end

	for column_xi = 1:J 
		if column_xi == row 
			values[1 + 4*K + column_xi + (row-1)*num_var] = (1/N)*sum{ tau[n,row]*(1-tau[n,row]) , n=1:N }
		else 
			values[1 + 4*K + column_xi + (row-1)*num_var] = (1/N)*sum{ tau[n,row]*tau[n,column_xi] , n=1:N }
	end

	for l = 1:L
		values[row + num_var - L + l + (row-1)*num_var] = 0 # first cluster of zeros
	end
end

# J*num_var values defined above. Now we start the "second line" of the jacobian, 0 -Z' eye(g)

for row2 = 1:L
	for col = 1:4*K+1
		values[col + J*num_var + (row2-1)*num_var] = 0 #second cluster of zeros
	end

	for col2 = 1:J 
		values[col2 + 4*K + 1 + J*num_var + (row2-1)*num_var] = -iv[row2, col2] # -Z'
	end 

	for col3 = 1:L 
		if col3 = row2
			values[col3 + num_var - L + J*num_var + (row2-1)*num_var] = g[col3] 
		else 																	# eye(g)
			values[col3 + num_var - L + J*num_var + (row2-1)*num_var] = 0
	end
end 


#################################
# "Put it all together version" #
#################################

for row = 1:J # each row is the derivative of s_j wrt (theta,xi);

	values[1 + (row-1)*num_var] = (1/N)*sum{ tau[n,row]*(p[row] - sum{ tau[n,j]*p[j] , j=1:J } ) , n=1:N} #derivatives of s_j w.r.t. alpha
	values[ 2 + (row-1)*num_var] = 0 # derivative w.r.t. beta_0
	values[2 + K + (row-1)*num_var] = (1/N)*sum{ tau[n,row]*( p[row] - sum{ tau[n,j]*p[j] , j=1:J } )*inc[n] , n=1:N } #piInc(1) => section that is multiplied by p
	values[2 + 2*K + (row-1)*num_var] = (1/N)*sum{ tau[n,row]*( p[row] - sum{ tau[n,j]*p[j] , j=1:J } )*age[n] , n=1:N } #piAge(1) => section that is multiplied by p
	values[2 + 3*K + (row-1)*num_var] = (1/N)*sum{ tau[n,row]*( p[row] - sum{ tau[n,j]*p[j] , j=1:J } )*v[n,1] , n=1:N } # w.r.t sigma_alpha


	for small_k = 2:K 
		values[1 + small_k + (row-1)*num_var] = (1/N)*sum{ tau[n, row]*( x[row, small_k] - sum{ tau[n,j]*x[j,small_k] , j=1:J} , n=1:N) } #w.r.t to beta_2, beta_3, beta_4
		values[1 + K + small_k + (row-1)*num_var] = (1/N)*sum{ tau[n, row]*(x[row,small_k] - sum{ tau[n,j]*x[j,small_k])*inc[n] , j=1:J } , n=1:N ) }
		#piInc(2), piInc(3) and piInc(4), the section multiplied by x
		values[1 + 2*K + small_k + (row-1)*num_var] = (1/N)*sum{ tau[n, row]*(x[row,small_k] - sum{ tau[n,j]*x[j,small_k])*age[n] , j=1:J } , n=1:N ) }
		#piAge(2), piAge(3) and piAge(4), the section multiplied by x
		values[1 + 3*K + small_k + (row-1)*num_var] = (1/N)*sum{ tau[n, row]*(x[row,small_k] - sum{ tau[n,j]*x[j,small_k])*v[n,small_k] , j=1:J } , n=1:N ) } 
		#sigma_beta_2, sigma_beta_3, sigma_beta_4
	end

	for column_xi = 1:J 
		if column_xi == row 
			values[1 + 4*K + column_xi + (row-1)*num_var] = (1/N)*sum{ tau[n,row]*(1-tau[n,row]) , n=1:N }
		else 
			values[1 + 4*K + column_xi + (row-1)*num_var] = (1/N)*sum{ tau[n,row]*tau[n,column_xi] , n=1:N }
	end

	for l = 1:L
		values[row + num_var - L + l] = 0 # first cluster of zeros
	end

	# J*num_var values defined above. Now we start the "second line" of the jacobian, 0 -Z' eye(g)

for row2 = 1:L
	for col = 1:4*K+1
		values[col + J*num_var + (row2-1)*num_var] = 0 #second cluster of zeros
	end

	for col2 = 1:J 
		values[col2 + 4*K + 1 + J*num_var + (row2-1)*num_var] = -iv[row2, col2] # -Z'
	end 

	for col3 = 1:L 
		if col3 = row2
			values[col3 + num_var - L + J*num_var + (row2-1)*num_var] = g[col3] 
		else 																	# eye(g)
			values[col3 + num_var - L + J*num_var + (row2-1)*num_var] = 0
	end
end