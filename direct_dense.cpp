// Direct exact evaluation of the specified Y_n or C_n, without a transfer.
// g++ -O3 -march=native -std=c++17 direct_dense.cpp -o direct_dense
// ./direct_dense 5 Y right
// No exterior-algebra formula, identity at a special arity, symmetry quotient,
// predicted coefficient, or multiplier is used in the recursion.
// Dense indexing is injective on valid subsets and merely saves storage.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#ifdef __BMI2__
#include <immintrin.h>
#endif
using Mask=uint64_t;
using Exp=std::array<int,5>;
struct Field{Exp alpha{};int derivative=0;};
static Mask extract(Mask s,Mask positions){
#ifdef __BMI2__
 return _pext_u64(s,positions);
#else
 Mask out=0,b=1;while(positions){Mask p=positions&(-positions);if(s&p)out|=b;positions^=p;b<<=1;}return out;
#endif
}
struct Solver{
 int n,m,nc,nq;bool right,isY,compact;
 std::vector<Field> fields;std::vector<Exp> delta;
 Mask diagonal=0,nondiagonal=0,full=0;
 std::vector<uint64_t> valid;
 std::vector<uint32_t> prefix;
 int32_t *memo=nullptr;
 static constexpr int32_t UNSET=std::numeric_limits<int32_t>::min();
 static constexpr int32_t BIG=UNSET+1;
 size_t cells=0,zero_words=0;
 uint64_t valid_bases=0,evaluations=0,zero_states=0,nonzero_states=0,
 zero_hits=0,nonzero_hits=0,max_abs=0,index_tests=0;
 int max_factor=0;
 std::unordered_map<Mask,int64_t> large;
 std::chrono::steady_clock::time_point start=std::chrono::steady_clock::now();
 double elapsed()const{return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
 Solver(int rank,bool y,bool r):n(rank),right(r),isY(y){
  nc=n-(isY?0:1);nq=n-(isY?1:2);
  for(int i=0;i<nc;++i){Field f;f.derivative=i;fields.push_back(f);}
  for(int a=0;a<n;++a)for(int b=0;b<n;++b){Field f;f.alpha[a]=1;f.derivative=b;if(a==b)diagonal|=Mask(1)<<fields.size();fields.push_back(f);}
  for(int i=0;i<nq;++i){Field f;f.alpha[i]=2;f.derivative=i;fields.push_back(f);}
  m=fields.size();full=(Mask(1)<<m)-1;nondiagonal=full^diagonal;
  for(auto f:fields){Exp d=f.alpha;--d[f.derivative];delta.push_back(d);}
  compact=m>23;
  if(compact)build_index();
  else {zero_words=((size_t(1)<<m)+63)/64;}
  cells=compact?(valid_bases<<n):(size_t(1)<<m);
  memo=static_cast<int32_t*>(malloc(cells*sizeof(int32_t)));
  if(!memo)throw std::runtime_error("allocation failure");
  std::fill(memo,memo+cells,UNSET);
  std::cerr<<"rank="<<n<<" family="<<(isY?'Y':'C')<<" arity="<<m<<" direction="<<(right?"right":"left")
   <<" valid_bases="<<valid_bases<<" memo_bytes="<<cells*4
   <<" index_bytes="<<valid.size()*8+prefix.size()*4<<"\n";
 }
 ~Solver(){free(memo);}
 static uint64_t hash(uint64_t x){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;x^=x>>31;return x;}
 void build_index(){
  // Remove diagonal linear fields only from the zero-bitset INDEX.
  // They remain distinct inputs in every recurrence and every subset key.
  // Their exponent shifts are zero, so validity is independent of their bits.
  const int ne=n*(n-1),basebits=nc+ne+nq;
  size_t basecount=size_t(1)<<basebits;
  valid.assign((basecount+63)/64,0);
  using Small=std::array<int8_t,5>;
  std::vector<Small> ds;
  for(int a=0;a<n;++a)for(int b=0;b<n;++b)if(a!=b){Small d{};d[a]=1;d[b]=-1;ds.push_back(d);}
  std::vector<Small> bal(size_t(1)<<ne);
  for(uint32_t a=1;a<bal.size();++a){unsigned i=__builtin_ctz(a);auto prev=bal[a&(a-1)];for(int j=0;j<n;++j)bal[a][j]=prev[j]+ds[i][j];}
  uint64_t pattern[32]{};
  for(unsigned constraint=0;constraint<(1u<<nc);++constraint)
   for(unsigned c=0;c<(1u<<nc);++c)
    if(right?((c&constraint)==constraint):((c&constraint)==0))pattern[constraint]|=uint64_t(1)<<c;
  for(unsigned q=0;q<(1u<<nq);++q)for(uint32_t a=0;a<bal.size();++a){
   bool ok=true;unsigned constraint=0;
   for(int i=0;i<n;++i){
    int qi=i<nq?((q>>i)&1):0;
    int v=bal[a][i]+qi+(right?0:(i==n-(isY?1:2)));
    if(right){ // mu_i=c_i-q_i-bal_i >=0
     if(v>(i<nc?1:0)){ok=false;break;}
     if(i<nc&&v==1)constraint|=1u<<i;
    }else{ // nu_i=lambda_i+q_i+bal_i-c_i >=0
     if(v<0){ok=false;break;}
     if(i<nc&&v==0)constraint|=1u<<i;
    }
   }
   if(ok){size_t first=(size_t(q)<<(ne+nc))|(size_t(a)<<nc);valid[first>>6]|=pattern[constraint]<<(first&63);}
  }
  prefix.resize(valid.size()+1);prefix[0]=0;
  for(size_t i=0;i<valid.size();++i)prefix[i+1]=prefix[i]+__builtin_popcountll(valid[i]);
  valid_bases=prefix.back();zero_words=((valid_bases<<n)+63)/64;
  // Independent checks of the validity table against the original field list.
  uint64_t seed=12983123;
  for(unsigned test=0;test<100000;++test){
   seed=hash(seed+0x9e3779b97f4a7c15ULL);Mask S=seed&full;
   Exp mu{};if(!right)mu[n-(isY?1:2)]=1;
   for(Mask bits=S;bits;bits&=bits-1){int i=__builtin_ctzll(bits);for(int j=0;j<n;++j)mu[j]+=right?-delta[i][j]:delta[i][j];}
   bool good=true;for(int j=0;j<n;++j)if(mu[j]<0)good=false;
   Mask base=extract(S,nondiagonal);bool indexed=(valid[base>>6]>>(base&63))&1;
   if(good!=indexed)throw std::runtime_error("validity index check failed");
   ++index_tests;
  }
 }
 uint64_t zero_index(Mask S)const{
  if(!compact)return S;
  Mask base=extract(S,nondiagonal);unsigned bit=base&63;size_t word=base>>6;
  if(((valid[word]>>bit)&1)==0)throw std::runtime_error("recursion reached invalid index");
  uint64_t rank=prefix[word]+__builtin_popcountll(valid[word]&((uint64_t(1)<<bit)-1));
  return (rank<<n)|extract(S,diagonal);
 }
 int64_t rec(Mask S,Exp mu){
  if(S==0)return 1;
  uint64_t zi=zero_index(S);
  int32_t saved=memo[zi];
  if(saved!=UNSET){
   if(saved==0)++zero_hits;else ++nonzero_hits;
   if(saved==BIG){auto it=large.find(S);if(it==large.end())throw std::runtime_error("large-value cache mismatch");return it->second;}
   return saved;
  }
  ++evaluations;
  if(evaluations%10000000==0)std::cerr<<"states="<<evaluations<<" zero="<<zero_states<<" nonzero="<<nonzero_states<<" seconds="<<elapsed()<<"\n";
  __int128 total=0;
  int sign=right&&(__builtin_popcountll(S)%2==0)?-1:1;
  for(Mask bits=S;bits;){
   Mask bit=bits&(-bits);bits^=bit;int i=__builtin_ctzll(bit);Exp nu=mu;bool good=true;
   for(int j=0;j<n;++j){nu[j]+=right?delta[i][j]:-delta[i][j];if(nu[j]<0)good=false;}
   int factor=right?mu[fields[i].derivative]:nu[fields[i].derivative];
   if(good&&factor>0){
    if(factor>(isY?n:n-1))throw std::runtime_error("derivative factor bound violated");
    max_factor=std::max(max_factor,factor);
    total+=__int128(sign)*factor*rec(S^bit,nu);
   }
   sign=-sign;
  }
  if(total<std::numeric_limits<int64_t>::min()||total>std::numeric_limits<int64_t>::max())throw std::overflow_error("exact value exceeds int64_t");
  int64_t value=int64_t(total);uint64_t mag=total<0?uint64_t(-total):uint64_t(total);max_abs=std::max(max_abs,mag);
  if(value==0)++zero_states;else ++nonzero_states;
  if(value>BIG && value<=std::numeric_limits<int32_t>::max())memo[zi]=int32_t(value);
  else{memo[zi]=BIG;large.emplace(S,value);}
  return value;
 }
 void run(){
  Exp lambda{};lambda[n-(isY?1:2)]=1;Exp out=lambda;
  for(auto d:delta)for(int j=0;j<n;++j)out[j]+=d[j];
  for(int j=0;j<n;++j)if(out[j]!=0)throw std::runtime_error("nonconstant final exponent");
  int64_t value=rec(full,right?lambda:out);
  std::cout<<"{\"rank\":"<<n<<",\"family\":\""<<(isY?'Y':'C')<<"\",\"arity\":"<<m
   <<",\"direction\":\""<<(right?"right":"left")<<"\",\"value\":"<<value
   <<",\"memo_states\":"<<evaluations<<",\"zero_states\":"<<zero_states<<",\"nonzero_states\":"<<nonzero_states
   <<",\"zero_hits\":"<<zero_hits<<",\"nonzero_hits\":"<<nonzero_hits<<",\"large_values\":"<<large.size()
   <<",\"max_abs_value\":"<<max_abs<<",\"max_factor\":"<<max_factor<<",\"index_tests\":"<<index_tests
   <<",\"seconds\":"<<elapsed()<<"}\n";
 }
};
int main(int argc,char**argv){try{
 int n=argc>1?std::stoi(argv[1]):5;std::string f=argc>2?argv[2]:"Y",d=argc>3?argv[3]:"right";
 if(n<2||n>5||(f!="Y"&&f!="C")||(d!="right"&&d!="left"))throw std::runtime_error("usage: direct_dense {2..5} {Y|C} {right|left}");
 Solver s(n,f=="Y",d=="right");s.run();return 0;
}catch(std::exception const&e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;}}
