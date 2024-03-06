function [DIAG_VEC, DIAG_LABELS, OFFDIAG_VEC, OFFDIAG_LABELS] = DIAG_OFFDIAG_SPLIT(S, p_vec)
	K = length(p_vec);
	DIAG_VEC = [];
	DIAG_LABELS = [];
	OFFDIAG_VEC = [];
	OFFDIAG_LABELS = [];
	for k = 1 : K
		for kp = 1 : K
			row_index = [(sum(p_vec(1 : (k - 1))) + 1) : sum(p_vec(1 : k))];
			col_index = [(sum(p_vec(1 : (kp - 1))) + 1) : sum(p_vec(1 : kp))];
			SUB_matrix = S(row_index, col_index);
			if (kp == k)
                idx = logical(triu(ones(size(SUB_matrix)), 1));
				v = SUB_matrix(idx);
				v_labels = transpose(repelem(strcat(num2str(kp), " diagonal block"), length(v)));
				idy = logical(tril(ones(size(SUB_matrix)), - 1));
				u = SUB_matrix(idy);
				u_labels = transpose(repelem(strcat(num2str(kp), " diagonal block"), length(u)));
				DIAG_VEC = [DIAG_VEC; v; u];
				DIAG_LABELS = [DIAG_LABELS; v_labels; u_labels];
            elseif (kp > k)
                w = SUB_matrix(:);
				w_labels = transpose(repelem(strcat(num2str(k), ", ", num2str(kp), " block"), length(w)));
				OFFDIAG_VEC = [OFFDIAG_VEC; w];
				OFFDIAG_LABELS = [OFFDIAG_LABELS; w_labels];
			else 
				OFFDIAG_VEC = [OFFDIAG_VEC];
				OFFDIAG_LABELS = [OFFDIAG_LABELS];
			end
		end
	end
end