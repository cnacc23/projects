import numpy as np

#to create a matrix that has all scores 
def max_scores(S, T, gap_pen, match, miss):

    matrix= [[[0, ''] for i in range(len(S)+1)] for j in range(len(T)+1)]
    
    #initialize first row and column scores 
    for i in range(0,len(T)):
        matrix[0][i][0]= i * gap_pen
    
    for i in range(0,len(S)):
        matrix[i][0][0]= i * gap_pen
    
   

    for i in range(1, len(S)):
        for j in range(1, len(T)):
           
            diag= matrix[i-1][j-1][0] 
               
            #diagonal score needs to compare bases to determine aln score 
            if S[i-1] ==  T[j-1]: 
                diag += match
            else:
                diag += miss 
              
            #up and side scores 
            up= matrix[i-1][j][0] + gap_pen
            left= matrix[i][j-1][0] + gap_pen
           
            #take max of three scores and assign it to cell[i, j]
            matrix[i][j][0]= max(diag, up, left)
            if max(diag, up, left) == diag:
                 matrix[i][j][1]= "d"
            elif max(diag, up, left) == up:
                 matrix[i][j][1]= "u"
            else:
                matrix[i][j][1]= "l"

    
    return matrix


#find optimal alignment of S and T 
def needleman_wunsch(S, T, matrix):

    #will give optimal alignment 
    opt_aln = [['' for i in range(len(S))] for j in range(len(T))]


    #start traceback at bottom corner of matrix 
    row= len(S)-1
    col= len(T)-1

    #current cell 
    curr= matrix[row][col][1]

    while row >= 0 and col >= 0:

        # if score comes from diagonal align two characters and move to diagonal cell
        if curr == "d":
            opt_aln[0][row]= S[row]
            opt_aln[1][col]= T[col]
            row -= 1
            col -= 1
    
        # if score comes from above, align char from S1 with - and move to cell above
        elif curr == "u":
            opt_aln[0][row]= S[row]
            opt_aln[1][col]= "-"
            row -= 1
    
        # if score comes from left, align char from S2 with - and move to left cell 
        elif curr == "l":
            opt_aln[0][row]= "-"
            opt_aln[1][col]= T[col]
            col -= 1
    
        # update curr 
        curr= matrix[row][col][1]

    for row in opt_aln:
        print(row)
         

"""
def score_assign(S, T, match, mismatch, cs, cn):

    #lengths of sequences 
    len_s= len(S)
    
   # align lists into a matrix 
    A= [S, T]

    # counter variable to keep track of optimal score 
    opt_score= 0

    # iterate through each sequence and compare same elements 
    for i in range(len_s):
    
            # bases match 
        if (A[0][i-1] == A[1][i-1]) and (A[0][i-1] != "-" and A[0][i-1] != "-"):
                opt_score += match 

            # bases don't match 
        elif (A[0][i-1] != A[1][i-1]) and (A[0][i-1] != "-" and A[1][i-1] != "-"):
                opt_score += mismatch 
                    
            # there is a gap in one of sequences         
        elif A[0][i-1] == "-":
            if A[1][i-1] == A[1][i-1]:
                opt_score += cs
            else:
                opt_score += cn
        elif A[1][i-1] == "-": 
            if A[0][i-1] == A[0][i-1]:
                opt_score += cs
            else:
                opt_score += cn

    #return opt score once calculations finished 
    print("Optimal score: ", opt_score)
"""

def main():
     
    S= ["A", "P", "P", "L", "E"]
    T= ["C", "A", "P", "E"]

    match= 1
    mismatch= -1
    cs= -1
    cn= -2

    
    score_matrix= max_scores(S, T, -1, match, mismatch)
    #needleman_wunsch(S, T, score_matrix)
    
    
if __name__ == "__main__":
     main()
