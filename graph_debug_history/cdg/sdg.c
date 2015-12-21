#include <stdio.h>
#include <stdlib.h>

#define EOS 0

typedef struct {
  int k; // length of node label
  int m; // number of nodes
  int p; // position of EOF
  int *last;
  int *W;
  int *F;
} sdg;

int u(int c)
{
  if (c >= 0) return c;
  else return -c;
}

int cc(int c)
{
  return (c == 0) ? '$' : c;
}

int C(sdg *G, int i)
{
int c;
  for (c = 0; c < 256; c++) {
    if (G->F[c] >= i) return c-1;
  }
  printf("??? C(%d)\n", i);
  return -1;
}

int *insert(int *S, int n, int i, int c)
{
  int *S2;
  int j;
  S2 = realloc(S, sizeof(int)*(n+2));
  for (j=n; j>=i; j--) S2[j+1] = S2[j];
  S2[i] = c;
  return S2;
}

int *delete(int *S, int n, int i)
{
  int *S2;
  int j;
  for (j=i; j<=n-1; j++) S[j] = S[j+1];
  S2 = realloc(S, sizeof(int)*(n));
  return S2;
}

int rank(int *S, int n, int i, int c)
{
int r, j;
  if (i <= 0) return 0;
  if (i > n) {
//    printf("rank: i = %d n = %d\n", i, n);
//    exit(1);
    i = n;
  }
  r = 0;
  for (j=1; j<=i; j++) {
    if (S[j] == c) r++;
  }
  return r;
}

int select2(int *S, int n, int i, int c)
{
int r, s;
  if (i <= 0) return 0;
  r = 0;
  for (s=1; s<=n; s++) {
    if (S[s] == c) r++;
    if (r == i) break;
  }
  return s;
}

int pred(int *S, int n, int i, int c)
{
  return select2(S, n, rank(S, n, i, c), c);
}

int succ(int *S, int n, int i, int c)
{
  return select2(S, n, rank(S, n, i-1, c)+1, c);
}

int range(sdg *G, int i, int *s, int *t)
{
  *s = pred(G->last, G->m, i-1, 1)+1;
  *t = succ(G->last, G->m, i, 1);
  return 0;
}

int bwd(sdg *G, int i)
{
  int c,r,j;
  c = C(G, i);
  r = rank(G->last, G->m, i, 1) - rank(G->last, G->m, G->F[c], 1);
  j = select2(G->W, G->m, r, c);
  return j;
}

int fwd(sdg *G, int j)
{
  int c,r,i;
  c = G->W[j];
  r = rank(G->W, G->m, j, c);
  i = select2(G->last, G->m, rank(G->last, G->m, G->F[c], 1) + r, 1);
  return i;
}

int cdeg(sdg *G, int v)
{
  if (G->last[v] != 1) {
    printf("cdeg: last[%d] = %d\n", v, G->last[v]);
    exit(1);
  }
  return v - pred(G->last, G->m, v-1, 1);
}

int child(sdg *G, int v, int c)
{
  int s,t;
  int j, j1, j2;
  range(G, v, &s, &t);
  j1 = pred(G->W, G->m, t, c);
  if (j1 >= s) {
    j = j1;
  } else {
    j2 = pred(G->W, G->m, t, -c);
    if (j2 >= s) {
      j = j1;
    } else {
      return -1;
    }
  }
  return fwd(G, j);
}

int pdeg(sdg *G, int v)
{
  int d,x,y,r;
  if (v == 1) return 0;
  d = C(G, v);
  x = bwd(G, v);
  y = succ(G->W, G->m, x+1, d);
  r = rank(G->W, G->m, y-1, -d) - rank(G->W, G->m, x+1, -d) + 1;
//  printf("pdeg d = %c x = %d y = %d\n", cc(d), x, y);
  return r;
}

int head(sdg *G, int v)
{
  int i, c;
//  printf("head %d ", v);
  if (v == 1) return EOS;
  for (i=0; i < G->k-1; i++) {
    v = succ(G->last, G->m, v, 1);
    v = bwd(G, v);
    if (v == 1) return EOS;
//    printf("%d ", v);
  }
  c = C(G, v);
//  printf("c = %c ", cc(c));
  return c;
}


int parent(sdg *G, int v, int c)
{
  int d,x,y;
  if (v == 1) return -1;
  d = C(G, v);
  x = bwd(G, v);
  if (head(G, x) == c) return succ(G->last, G->m, x, 1);
  y = succ(G->W, G->m, x+1, d);
  while (1) {
    x = succ(G->W, G->m, x, -d);
    if (x >= y) break;
    if (head(G, x) == c) return succ(G->last, G->m, x, 1);
  }
  return -1;
}

void print_sdg(sdg *G)
{
  int i,j,d;
  int s,t;

  printf("k = %d m = %d\n", G->k, G->m);
  for (i=1; i <= G->m; i++) {
    printf("%3d: %d  ", i, G->last[i]);
    printf("|%c ", cc(head(G, i)));
    printf("%c| ", cc(C(G, i)));
    if (G->W[i] == EOS) {
      printf("$ ");
    } else if (G->W[i] > 0) {
      printf("%c ", G->W[i]);
    } else if (G->W[i] < 0) {
      printf("%c-", -G->W[i]);
    }
    range(G, i, &s, &t);
    printf("[%2d,%2d]", s, t);
    if (G->last[i] == 1) printf("(%2d)", bwd(G, i));
    else printf("(  )");
    printf("(%2d)", fwd(G, i));
//    printf("%3d", fwd(G, i));
    if (G->last[i] == 1) {
      d = cdeg(G, i);
      printf(" c:%d", d);
      for (j=d-1; j>=0; j--) {
        printf(" (%c:%d)", cc(u(G->W[i-j])), child(G, i, u(G->W[i-j])));
      }
    }
    if (G->last[i] == 1) {
      d = pdeg(G, i);
      printf(" p:%d", d);
    }
    printf("\n");
  }
  printf("F\n");
  for (i=0; i<256-1; i++) {
    if (G->F[i] < G->F[i+1]) {
      printf("%c %d\n", cc(i), G->F[i]);
    }
  }
  printf("\n");
}


void update_sdg(sdg *G, int c)
{
  int p,pp;
  int s,t;
  int i,j,j2,j3;
  int x,y,x2;
  int f1, f2, f3;
  int c2;
  p = G->p;
  if (G->W[p] != EOS) {
    printf("updat_sdg: p = %d W = %c\n", p, G->W[p]);
    exit(1);
  }
  f1 = 1;
  f2 = 1;
  f3 = 1;
  range(G, p, &s, &t);
  pp = succ(G->W, G->m, p, c);
  if (pp <= t) {
    f1 = 0;   // c aleady exists in R(p)
  } else {
    pp = succ(G->W, G->m, p, -c);
    if (pp <= t) { // c- aleady exists in R(p)
      f1 = 0;
    }
  }
  if (f1 == 0) { // c or c- aleady exists in R(p)
    c2 = C(G, p);
//    G->W = delete(G->W, G->m, p);
//    G->last = delete(G->last, G->m, p);
//    G->m--;
//    for (i=c2+1; i<256; i++) G->F[i]--;

    x = fwd(G, pp);

    G->W = delete(G->W, G->m, p);
    G->last = delete(G->last, G->m, p);

    if (p < x) x--;

    G->W = insert(G->W, G->m-1, x, EOS);
    G->last = insert(G->last, G->m-1, x, 0);
    G->p = x;
//    G->m++;
    for (i=c2+1; i<256; i++) G->F[i]--;
    for (i=c+1; i<256; i++) G->F[i]++;

  } else { // c does not exist
    j = pred(G->W, G->m, p-1, c);
    if (j > 0) {
      x = fwd(G, j);
      // $ will be inserted at W[x+1]
      if (G->last[x] == 0) { // Node[x] == Node[x+1]
        f2 = 0; // node label of $ is the same as W[x]
      } else {
        j2 = succ(G->W, G->m, p+1, c);
        j3 = succ(G->W, G->m, p+1, -c);
        if (j3 < j2) { // Node'[j3] == Node'[j]
          f2 = 0; // node label of $ is the same as W[x]
        } else { // check if Node[j] == Node[p]
          int i1, i2;
          i1 = j;  i2 = p;
          f2 = 0;
          for (i=0; i < G->k-1; i++) {
            if (C(G, i1) != C(G, i2)) {
              f2 = 1;
              break;
            }
            i1 = bwd(G, i1);
            i2 = bwd(G, i2);
          }
          if (j2 <= G->m) {
            i1 = j2;  i2 = p;
            f3 = 0;
            y = j2;
            for (i=0; i < G->k-1; i++) {
              if (C(G, i1) != C(G, i2)) {
                f3 = 1;
                break;
              }
              i1 = bwd(G, i1);
              i2 = bwd(G, i2);
            }
          }
        }
      }

      if (f2 == 0) {
        G->W[p] = -c; // $ -> c
      } else {
        G->W[p] = c; // $ -> c
      }
      if (f3 == 0) {
        G->W[y] = -c;
      }

      if (f2 == 0 || f3 == 0) {
        G->W = insert(G->W, G->m, x, EOS);
        G->last = insert(G->last, G->m, x, 0);
        G->p = x;
      } else {
        G->W = insert(G->W, G->m, x+1, EOS);
        G->last = insert(G->last, G->m, x+1, 1);
        G->p = x+1;
      }
      G->m++;
      for (i=c+1; i<256; i++) G->F[i]++;
    } else {
      int i1, i2;
      x = G->F[c];
      // f2 = 1;
      j2 = succ(G->W, G->m, p+1, c);
      if (j2 <= G->m) {
        i1 = j2;  i2 = p;
        f3 = 0;
        y = j2;
        for (i=0; i < G->k-1; i++) {
          if (C(G, i1) != C(G, i2)) {
            f3 = 1;
            break;
          }
          i1 = bwd(G, i1);
          i2 = bwd(G, i2);
        }
      }
      if (f2 == 0) {
        G->W[p] = -c; // $ -> c
      } else {
        G->W[p] = c; // $ -> c
      }
      if (f3 == 0) {
        G->W[y] = -c;
      }

      if (f2 == 0 || f3 == 0) {
        G->W = insert(G->W, G->m, x+1, EOS);
        G->last = insert(G->last, G->m, x+1, 0);
        G->p = x+1;
      } else {
        G->W = insert(G->W, G->m, x+1, EOS);
        G->last = insert(G->last, G->m, x+1, 1);
        G->p = x+1;
      }
      G->m++;
      for (i=c+1; i<256; i++) G->F[i]++;
    }

  }
  

}

sdg *init_sdg(int k)
{
sdg *G;
int i, m;
  G = malloc(sizeof(sdg));
  if (G == NULL) {
    printf("init_sdg: malloc\n");
    exit(1);
  }
  G->k = k;
  G->m = m = 1;
  G->last = malloc(sizeof(*G->last) * (m+1));
  G->W = malloc(sizeof(*G->W) * (m+1));
  G->F = malloc(sizeof(*G->F) * 256);
  if (G->last == NULL || G->W == NULL || G->F == NULL) {
    printf("init_sdg: malloc\n");
    exit(1);
  }
  G->last[m] = 1;
  G->W[m] = EOS;
  G->p = m;

  G->F[0] = 0;
  for (i=1; i<256; i++) G->F[i] = 1;
  return G;
}

sdg *init_sdg2(int k, int m, int *last, int *W, int *W2)
{
sdg *G;
int i,j,r;
  G = malloc(sizeof(sdg));
  if (G == NULL) {
    printf("init_sdg: malloc\n");
    exit(1);
  }
  G->k = k;
  G->m = m;
  G->last = malloc(sizeof(*G->last) * (m+1));
  G->W = malloc(sizeof(*G->W) * (m+1));
  G->F = malloc(sizeof(*G->F) * 256);
  if (G->last == NULL || G->W == NULL || G->F == NULL) {
    printf("init_sdg2: malloc\n");
    exit(1);
  }
  for (i=1; i<=m; i++) {
    G->last[i] = last[i];
    G->W[i] = W[i];
    if (W[i] == EOS) G->p = i;
  }
  for (i=0; i<256; i++) G->F[i] = 0;
  for (i=1; i<=m; i++) {
    G->F[W2[i]]++;
  }
  r = 0;
  for (i=0; i<256; i++) {
    j = G->F[i];
    G->F[i] = r;
    r += j;
  }

  return G;
}

sdg *test0()
{
sdg *G;
  int last[13] = {-1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1};
  int W[14] = {0, 'T', EOS, 'C', 'C', 'G', -'G', 'G', 'A', 'T', -'A', 'A', 'C'};
  int W2[14] = {0, EOS, 'A', 'A', 'A', 'C', 'C', 'C', 'G', 'G', 'G', 'T', 'T'};
  G = init_sdg2(3, 12, last, W, W2);
  print_sdg(G);
  return G;
}

sdg *test()
{
sdg *G;
  int last[14] = {-1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1};
  int W[14] = {0, 'T', 'C', 'C', 'G', 'T', -'G', 'G', 'A', 'T', -'A', 'A', EOS, 'C'};
  int W2[14] = {0, EOS, 'A', 'A', 'C', 'C', 'C', 'C', 'G', 'G', 'G', 'T', 'T', 'T'};
  G = init_sdg2(3, 13, last, W, W2);
  print_sdg(G);
  return G;
}


int main(int argc, char *argv[])
{
sdg *G;
//int T[13] = {'T', 'A', 'C', 'G', 'A', 'C', 'G', 'T', 'C', 'G', 'A', 'C', 'T'};
int T[7] = {'A', 'A', 'A', 'A', 'A', 'A', 'A'};
int i;

#if 0
  G = test();

  update_sdg(G, 'C');
  print_sdg(G);

  update_sdg(G, 'T');
  print_sdg(G);
#endif

  G = init_sdg(3);
  print_sdg(G);
  
  for (i=0; i<13; i++) {
    update_sdg(G, T[i]);
    print_sdg(G);
  }

}
