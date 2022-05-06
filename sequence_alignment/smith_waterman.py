import numpy as np

def smith_waterman(seq1, seq2, match, mismatch, gapopen, gapextend):
    max_score = 0
    M = {}
    S = np.zeros([len(seq1) + 1, len(seq2) + 1])
    
    for i in range(0,len(seq1)+1):
        for j in range(0,len(seq2)+1):
            
            M['[' + str(i) + ', ' + str(j) + ']'] = ""
            
            if i == 0 or j == 0:
                pass
            else:
                if seq1[i-1] == seq2[j-1]:
                    s = S[i-1,j-1] + match
                else:
                    s = S[i-1,j-1] - mismatch
                d = S[i-1,j-1] + s
                
                if M['[' + str(i-1) + ', ' + str(j) + ']'] != "u":
                    p = gapopen + gapextend
                else:
                    p = gapextend            
                u = S[i-1,j] - p
                
                if M['[' + str(i) + ', ' + str(j-1) + ']'] != "l":
                    p = gapopen + gapextend            
                else:
                    p = gapextend                
                l = S[i,j-1] - p
                
                S[i,j] = max(s,u,l,0)
                
                
                if s == S[i,j]:
                    M['[' + str(i) + ', ' + str(j) + ']'] = "d"
                elif u == S[i,j]:
                    M['[' + str(i) + ', ' + str(j) + ']'] = "u"
                elif l == S[i,j]:
                    M['[' + str(i) + ', ' + str(j) + ']'] = "l"
                else:
                    M['[' + str(i) + ', ' + str(j) + ']'] = ""
                #print(i,j,S[i,j],s,u,l,M['[' + str(i) + ', ' + str(j) + ']'])
                     
                    
            if S[i,j] >= max_score:
                max_score = int(S[i,j])
                max_i = i
                max_j = j
	print(S)
    ## trace back
    alnseq1 = ""
    alnseq2 = ""

    m = max_i #index the base of seq1 
    n = max_j #index the base of seq2
    
    i = 0 #index the base of alnseq1
    j = 0 #index the base of alnseq2
    
    
    while 1:
        m = max_i - i
        n = max_j - j
        
        if M['[' + str(m) + ', ' + str(n) + ']'] == "d":
            alnseq1 = seq1[m-1] + alnseq1
            alnseq2 = seq2[n-1] + alnseq2
            i = i+1
            j = j+1
        if M['[' + str(m) + ', ' + str(n) + ']'] == "u":
            alnseq1 = seq1[m-1] + alnseq1
            alnseq2 = "-" + alnseq2
            i = i+1
        if M['[' + str(m) + ', ' + str(n) + ']'] == "l":
            alnseq1 = "-" + alnseq1
            alnseq2 = seq2[n-1] + alnseq2
            j = j+1
        #print(m, n, S[m,n],M['[' + str(m) + ', ' + str(n) + ']'],alnseq1,alnseq2)
        
        if M['[' + str(m) + ', ' + str(n) + ']'] == "":
            break

    return max_score, alnseq1, alnseq2
    
