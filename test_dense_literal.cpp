// Checks the direct evaluator against the definition as a permutation sum.
#define main direct_cli_main
#include "direct_dense.cpp"
#undef main
#include <numeric>
#include <random>

int64_t literal(const std::vector<Field>& fields, Exp lambda, int n){
 std::vector<int> p(fields.size());std::iota(p.begin(),p.end(),0);
 int64_t sum=0;
 do{
  int sign=1;for(size_t i=0;i<p.size();++i)for(size_t j=i+1;j<p.size();++j)if(p[i]>p[j])sign=-sign;
  Exp mu=lambda;int64_t c=sign;
  for(auto it=p.rbegin();it!=p.rend();++it){
   const auto&f=fields[*it];c*=mu[f.derivative];if(c==0)break;
   for(int j=0;j<n;++j)mu[j]+=f.alpha[j];
   --mu[f.derivative];
  }
  sum+=c;
 }while(std::next_permutation(p.begin(),p.end()));return sum;
}
int64_t evaluate(const std::vector<Field>& f, Exp lambda,bool right){
 Solver s(3,true,right); // allocation sufficient for all test masks
 s.fields=f;s.m=f.size();s.full=(Mask(1)<<s.m)-1;s.delta.clear();
 Exp out=lambda;
 for(auto x:f){Exp d=x.alpha;--d[x.derivative];s.delta.push_back(d);for(int j=0;j<3;++j)out[j]+=d[j];}
 if(!right)for(int j=0;j<3;++j)if(out[j]<0)return 0;
 return s.rec((Mask(1)<<s.m)-1,right?lambda:out);
}
int main(){try{
 std::mt19937 rng(194931);int count=0,nonzero=0,repeat_checks=0;
 for(int t=0;t<160;++t){
  int m=2+rng()%7;std::vector<Field> fs;
  int nq=0;
  for(int i=0;i<m;++i){Field f;f.derivative=rng()%3;int degree=rng()%3;if(degree==2&&nq==2)degree=1;if(degree==2)++nq;for(int j=0;j<degree;++j)++f.alpha[rng()%3];fs.push_back(f);}
  Exp lam{};lam[rng()%3]=1;
  int64_t expected=literal(fs,lam,3);
  for(bool right:{false,true}){
   int64_t got=evaluate(fs,lam,right);
   if(got!=expected)throw std::runtime_error("literal mismatch");
   // Reconstruct the solver and repeat; this checks repeatability, not collisions.
   if(evaluate(fs,lam,right)!=expected)throw std::runtime_error("repeat check mismatch");
   ++repeat_checks;
  }
  ++count;nonzero+=expected!=0;
 }
 // Fixed nonzero convention test: s4(D1,D2,x1D1,x2D1)(x1^2)=-4.
 std::vector<Field> f(4);f[0].derivative=0;f[1].derivative=1;
 f[2].derivative=0;f[2].alpha[0]=1;f[3].derivative=0;f[3].alpha[1]=1;
 Exp lam{};lam[0]=2;
 if(literal(f,lam,3)!=-4||evaluate(f,lam,false)!=-4||evaluate(f,lam,true)!=-4)throw std::runtime_error("s4 convention mismatch");
 std::cout<<"{\"literal_cases\":"<<count<<",\"nonzero_cases\":"<<nonzero
  <<",\"repeat_checks\":"<<repeat_checks<<",\"s4_convention\":-4,\"passed\":true}\n";
 return 0;
}catch(std::exception const&e){std::cerr<<e.what()<<"\n";return 1;}}
