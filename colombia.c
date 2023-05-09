#include <stdio.h>
#include <stdlib.h> 


int* max_num_cafes(int n, int ratings[], int *num_visited_cafes) {
    int dp[n];  //max number of cafes that can be visited 
    int last_visited_cafe[n];


    //populate dp with 1s, represents every cafe being visited alone 
    for (int i = 0; i < n; i++) {
        dp[i] = 1;
        last_visited_cafe[i] = -1;
    }


    int max_dp = 1;
    int last_visited_index = 0;
    
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            if (ratings[i] >= ratings[j] && dp[i] < dp[j] + 1) {
                dp[i] = dp[j] + 1;
                last_visited_cafe[i] = j;
                if (dp[i] > max_dp) {
                    max_dp = dp[i];
                    last_visited_index = i;
                }
            }
        }
    }

    //keep track of all cafes visited 
    *num_visited_cafes = max_dp;
    int *visited = malloc(max_dp * sizeof(int));
    int index = max_dp - 1;
    int cafe_index = last_visited_index;
    while (index >= 0) {
        visited[index] = ratings[cafe_index];
        cafe_index = last_visited_cafe[cafe_index];
        index--;
    }
    return visited;
}



int enjoyment_score(int n, int visited_cafes [], int num_visited){
    int first= visited_cafes[0]; 
    int last= visited_cafes[num_visited-1];   


    return last-first; 
    
}


int main(){

int n;
scanf("%d", &n);

//scan in n integers 
int ratings[n]; 
for(int i=0; i<n; i++){
    scanf("%d", &ratings[i]); 
}

int num_visited_cafes;

//calculate number of cafes that can be visited
int *visited_cafes = max_num_cafes(n, ratings, &num_visited_cafes);

for(int i=0; i<n; i++){
    printf("%d ", visited_cafes[i]);
}

free(visited_cafes);

//calculate difference in enoyment scores 
int diff_enjoyment= enjoyment_score(n, visited_cafes, num_visited_cafes); 

printf("\n"); 
printf("%d %d\n", num_visited_cafes, diff_enjoyment);

return 0; 

}
