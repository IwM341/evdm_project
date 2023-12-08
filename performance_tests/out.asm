	.file	"sum_vecto_test.cpp"
 # GNU C++23 (x86_64-posix-seh-rev3, Built by MinGW-W64 project) version 11.2.0 (x86_64-w64-mingw32)
 #	compiled by GNU C version 11.2.0, GMP version 6.2.1, MPFR version 4.1.0, MPC version 1.2.1, isl version isl-0.24-GMP

 # GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
 # options passed: -march=alderlake -mmmx -mpopcnt -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -mavx2 -mno-sse4a -mno-fma4 -mno-xop -mfma -mbmi -mbmi2 -maes -mpclmul -mno-avx512vl -mno-avx512bw -mno-avx512dq -mno-avx512cd -mno-avx512er -mno-avx512pf -mno-avx512vbmi -mno-avx512ifma -mno-avx5124vnniw -mno-avx5124fmaps -mno-avx512vpopcntdq -mno-avx512vbmi2 -mgfni -mvpclmulqdq -mno-avx512vnni -mno-avx512bitalg -mno-avx512bf16 -mno-avx512vp2intersect -mno-3dnow -madx -mabm -mno-cldemote -mclflushopt -mclwb -mno-clzero -mcx16 -mno-enqcmd -mf16c -mfsgsbase -mfxsr -mno-hle -msahf -mno-lwp -mlzcnt -mmovbe -mmovdir64b -mmovdiri -mno-mwaitx -mno-pconfig -mno-pku -mno-prefetchwt1 -mprfchw -mptwrite -mrdpid -mrdrnd -mrdseed -mno-rtm -mserialize -mno-sgx -msha -mshstk -mno-tbm -mno-tsxldtrk -mvaes -mwaitpkg -mno-wbnoinvd -mxsave -mxsavec -mxsaveopt -mxsaves -mno-amx-tile -mno-amx-int8 -mno-amx-bf16 -mno-uintr -mhreset -mno-kl -mno-widekl -mavxvnni --param=l1-cache-size=48 --param=l1-cache-line-size=64 --param=l2-cache-size=24576 -mtune=alderlake -mavx512f -O3 -std=c++23
	.text
	.p2align 4
	.def	__tcf_1;	.scl	3;	.type	32;	.endef
	.seh_proc	__tcf_1
__tcf_1:
.LFB14360:
	.seh_endprologue
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/iostream:74:   static ios_base::Init __ioinit;
	leaq	_ZStL8__ioinit(%rip), %rcx	 #, tmp82
	jmp	_ZNSt8ios_base4InitD1Ev	 #
	.seh_endproc
	.p2align 4
	.globl	_Z8sum_simdPKfy
	.def	_Z8sum_simdPKfy;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z8sum_simdPKfy
_Z8sum_simdPKfy:
.LFB13472:
	subq	$120, %rsp	 #,
	.seh_stackalloc	120
	.seh_endprologue
 # sum_vecto_test.cpp:20: 	std::memcpy(tmp,poly,16*sizeof(float));
	vmovdqa	(%rcx), %xmm5	 # MEM <unsigned char[64]> [(char * {ref-all})_26], tmp176
 # sum_vecto_test.cpp:17: float __attribute__ ((noinline)) sum_simd(const float * array,const size_t N){
	leaq	63(%rsp), %r8	 #, tmp126
	andq	$-64, %r8	 #, tmp128
 # sum_vecto_test.cpp:20: 	std::memcpy(tmp,poly,16*sizeof(float));
	vmovdqa	%xmm5, (%r8)	 # tmp176, MEM <unsigned char[64]> [(char * {ref-all})&tmp]
	vmovdqa	16(%rcx), %xmm5	 # MEM <unsigned char[64]> [(char * {ref-all})_26], tmp177
	vmovdqa	32(%rcx), %xmm4	 # MEM <unsigned char[64]> [(char * {ref-all})_26], tmp178
	vmovdqa	%xmm5, 16(%r8)	 # tmp177, MEM <unsigned char[64]> [(char * {ref-all})&tmp]
	vmovdqa	48(%rcx), %xmm5	 # MEM <unsigned char[64]> [(char * {ref-all})_26], tmp179
 # sum_vecto_test.cpp:17: float __attribute__ ((noinline)) sum_simd(const float * array,const size_t N){
	movq	%rdx, %r9	 # N, tmp174
 # sum_vecto_test.cpp:20: 	std::memcpy(tmp,poly,16*sizeof(float));
	vmovdqa	%xmm4, 32(%r8)	 # tmp178, MEM <unsigned char[64]> [(char * {ref-all})&tmp]
	vmovdqa	%xmm5, 48(%r8)	 # tmp179, MEM <unsigned char[64]> [(char * {ref-all})&tmp]
 # sum_vecto_test.cpp:22: 	for(size_t i=1;i<N/16;++i){
	shrq	$4, %rdx	 #, _11
 # sum_vecto_test.cpp:22: 	for(size_t i=1;i<N/16;++i){
	cmpq	$31, %r9	 #, tmp174
	jbe	.L4	 #,
	salq	$6, %rdx	 #, tmp136
	vmovaps	(%r8), %zmm0	 # MEM <vector(16) float> [(float *)&tmp], vect_tmp_I_lsm.125
	leaq	64(%rcx), %rax	 #, ivtmp.131
	addq	%rcx, %rdx	 # array, _137
	.p2align 4,,10
	.p2align 3
.L5:
 # sum_vecto_test.cpp:24: 			tmp[j] += poly[i*16+j];
	vaddps	(%rax), %zmm0, %zmm0	 # MEM <const vector(16) float> [(const float *)_148], vect_tmp_I_lsm.125, vect_tmp_I_lsm.125
	addq	$64, %rax	 #, ivtmp.131
	cmpq	%rax, %rdx	 # ivtmp.131, _137
	jne	.L5	 #,
	vmovaps	%zmm0, (%r8)	 # vect_tmp_I_lsm.125, MEM <vector(16) float> [(float *)&tmp]
.L4:
 # sum_vecto_test.cpp:29: 		sum += tmp[j];
	vmovaps	(%r8), %zmm1	 # MEM <vector(16) float> [(float *)&tmp], vect__9.110
	vxorps	%xmm3, %xmm3, %xmm3	 # stmp_sum_21.111
	vaddss	%xmm1, %xmm3, %xmm3	 # tmp138, stmp_sum_21.111, stmp_sum_21.111
	vshufps	$85, %xmm1, %xmm1, %xmm4	 #, tmp138, tmp138, stmp_sum_21.111
	vshufps	$255, %xmm1, %xmm1, %xmm2	 #, tmp138, tmp138, stmp_sum_21.111
	vaddss	%xmm3, %xmm4, %xmm4	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vunpckhps	%xmm1, %xmm1, %xmm3	 # tmp138, tmp138, stmp_sum_21.111
	vextractf128	$0x1, %ymm1, %xmm0	 # tmp137, tmp147
	vaddss	%xmm4, %xmm3, %xmm3	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vextracti64x4	$0x1, %zmm1, %ymm1	 # vect__9.110, tmp154
	vaddss	%xmm3, %xmm2, %xmm2	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vshufps	$85, %xmm0, %xmm0, %xmm3	 #, tmp147, tmp147, stmp_sum_21.111
	vaddss	%xmm2, %xmm0, %xmm2	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vaddss	%xmm2, %xmm3, %xmm3	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vunpckhps	%xmm0, %xmm0, %xmm2	 # tmp147, tmp147, stmp_sum_21.111
	vshufps	$255, %xmm0, %xmm0, %xmm0	 #, tmp147, tmp147, stmp_sum_21.111
	vaddss	%xmm3, %xmm2, %xmm2	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vshufps	$85, %xmm1, %xmm1, %xmm3	 #, tmp155, tmp155, stmp_sum_21.111
	vaddss	%xmm2, %xmm0, %xmm0	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vshufps	$255, %xmm1, %xmm1, %xmm2	 #, tmp155, tmp155, stmp_sum_21.111
	vaddss	%xmm0, %xmm1, %xmm0	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vaddss	%xmm0, %xmm3, %xmm3	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vunpckhps	%xmm1, %xmm1, %xmm0	 # tmp155, tmp155, stmp_sum_21.111
	vextractf128	$0x1, %ymm1, %xmm1	 # tmp154, tmp163
	vaddss	%xmm3, %xmm0, %xmm0	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vaddss	%xmm0, %xmm2, %xmm2	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vaddss	%xmm2, %xmm1, %xmm0	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vshufps	$85, %xmm1, %xmm1, %xmm2	 #, tmp163, tmp163, stmp_sum_21.111
	vaddss	%xmm0, %xmm2, %xmm2	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vunpckhps	%xmm1, %xmm1, %xmm0	 # tmp163, tmp163, stmp_sum_21.111
 # sum_vecto_test.cpp:29: 		sum += tmp[j];
	vshufps	$255, %xmm1, %xmm1, %xmm1	 #, tmp163, tmp163, stmp_sum_21.111
	vaddss	%xmm2, %xmm0, %xmm0	 # stmp_sum_21.111, stmp_sum_21.111, stmp_sum_21.111
	vaddss	%xmm1, %xmm0, %xmm0	 # stmp_sum_21.111, stmp_sum_21.111, sum
	vzeroupper
 # sum_vecto_test.cpp:32: }
	addq	$120, %rsp	 #,
	ret	
	.seh_endproc
	.p2align 4
	.globl	_Z3sumPKfy
	.def	_Z3sumPKfy;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z3sumPKfy
_Z3sumPKfy:
.LFB13477:
	.seh_endprologue
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	testq	%rdx, %rdx	 # N
	je	.L15	 #,
	leaq	-1(%rdx), %rax	 #, tmp174
	cmpq	$14, %rax	 #, tmp174
	jbe	.L16	 #,
	movq	%rdx, %r8	 # N, bnd.139
	shrq	$4, %r8	 #, bnd.139
	salq	$6, %r8	 #, tmp176
	movq	%rcx, %rax	 # array, ivtmp.157
	addq	%rcx, %r8	 # array, _12
 # sum_vecto_test.cpp:35: 	float sum = 0;
	vxorps	%xmm0, %xmm0, %xmm0	 # <retval>
	.p2align 4,,10
	.p2align 3
.L10:
	vmovups	(%rax), %ymm1	 # MEM <const vector(16) float> [(const float *)_77], tmp178
	vmovups	(%rax), %zmm4	 # MEM <const vector(16) float> [(const float *)_77], tmp235
	vaddss	%xmm1, %xmm0, %xmm0	 # stmp_sum_10.145, <retval>, stmp_sum_10.145
	vshufps	$85, %xmm1, %xmm1, %xmm3	 #, tmp179, tmp179, stmp_sum_10.145
	vshufps	$255, %xmm1, %xmm1, %xmm2	 #, tmp179, tmp179, stmp_sum_10.145
	vaddss	%xmm3, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	vunpckhps	%xmm1, %xmm1, %xmm3	 # tmp179, tmp179, stmp_sum_10.145
	vextractf128	$0x1, %ymm1, %xmm1	 # tmp178, tmp187
	vaddss	%xmm3, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	addq	$64, %rax	 #, ivtmp.157
	vaddss	%xmm2, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	vshufps	$85, %xmm1, %xmm1, %xmm2	 #, tmp187, tmp187, stmp_sum_10.145
	vaddss	%xmm1, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	vaddss	%xmm2, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	vunpckhps	%xmm1, %xmm1, %xmm2	 # tmp187, tmp187, stmp_sum_10.145
	vshufps	$255, %xmm1, %xmm1, %xmm1	 #, tmp187, tmp187, stmp_sum_10.145
	vaddss	%xmm2, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	vaddss	%xmm1, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	vextracti64x4	$0x1, %zmm4, %ymm1	 # tmp235, tmp194
	vshufps	$85, %xmm1, %xmm1, %xmm3	 #, tmp195, tmp195, stmp_sum_10.145
	vaddss	%xmm1, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	vshufps	$255, %xmm1, %xmm1, %xmm2	 #, tmp195, tmp195, stmp_sum_10.145
	vaddss	%xmm3, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	vunpckhps	%xmm1, %xmm1, %xmm3	 # tmp195, tmp195, stmp_sum_10.145
	vextractf128	$0x1, %ymm1, %xmm1	 # tmp194, tmp203
	vaddss	%xmm3, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	vaddss	%xmm2, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	vshufps	$85, %xmm1, %xmm1, %xmm2	 #, tmp203, tmp203, stmp_sum_10.145
	vaddss	%xmm1, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	vaddss	%xmm2, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	vunpckhps	%xmm1, %xmm1, %xmm2	 # tmp203, tmp203, stmp_sum_10.145
 # sum_vecto_test.cpp:37: 		sum += array[i];
	vshufps	$255, %xmm1, %xmm1, %xmm1	 #, tmp203, tmp203, stmp_sum_10.145
	vaddss	%xmm2, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, stmp_sum_10.145
	vaddss	%xmm1, %xmm0, %xmm0	 # stmp_sum_10.145, stmp_sum_10.145, <retval>
	cmpq	%rax, %r8	 # ivtmp.157, _12
	jne	.L10	 #,
	movq	%rdx, %rax	 # N, tmp.149
	andq	$-16, %rax	 #, tmp.149
	testb	$15, %dl	 #, N
	je	.L24	 #,
.L9:
	movq	%rdx, %r8	 # N, niters.146
	subq	%rax, %r8	 # tmp.149, niters.146
	leaq	-1(%r8), %r9	 #, tmp212
	cmpq	$6, %r9	 #, tmp212
	jbe	.L13	 #,
 # sum_vecto_test.cpp:37: 		sum += array[i];
	vmovups	(%rcx,%rax,4), %ymm1	 # MEM <const vector(8) float> [(const float *)vectp_array.151_85], MEM <const vector(8) float> [(const float *)vectp_array.151_85]
	movq	%r8, %r9	 # niters.146, niters_vector_mult_vf.148
	vaddss	%xmm1, %xmm0, %xmm0	 # stmp_sum_7.153, <retval>, stmp_sum_7.153
	vshufps	$85, %xmm1, %xmm1, %xmm3	 #, tmp214, tmp214, stmp_sum_7.153
	vshufps	$255, %xmm1, %xmm1, %xmm2	 #, tmp214, tmp214, stmp_sum_7.153
	vaddss	%xmm3, %xmm0, %xmm0	 # stmp_sum_7.153, stmp_sum_7.153, stmp_sum_7.153
	vunpckhps	%xmm1, %xmm1, %xmm3	 # tmp214, tmp214, stmp_sum_7.153
	vextractf128	$0x1, %ymm1, %xmm1	 # MEM <const vector(8) float> [(const float *)vectp_array.151_85], tmp218
	vaddss	%xmm3, %xmm0, %xmm0	 # stmp_sum_7.153, stmp_sum_7.153, stmp_sum_7.153
	andq	$-8, %r9	 #, niters_vector_mult_vf.148
	addq	%r9, %rax	 # niters_vector_mult_vf.148, tmp.149
	vaddss	%xmm2, %xmm0, %xmm0	 # stmp_sum_7.153, stmp_sum_7.153, stmp_sum_7.153
	vshufps	$85, %xmm1, %xmm1, %xmm2	 #, tmp218, tmp218, stmp_sum_7.153
	vaddss	%xmm1, %xmm0, %xmm0	 # stmp_sum_7.153, stmp_sum_7.153, stmp_sum_7.153
	vaddss	%xmm2, %xmm0, %xmm0	 # stmp_sum_7.153, stmp_sum_7.153, stmp_sum_7.153
	vunpckhps	%xmm1, %xmm1, %xmm2	 # tmp218, tmp218, stmp_sum_7.153
 # sum_vecto_test.cpp:37: 		sum += array[i];
	vshufps	$255, %xmm1, %xmm1, %xmm1	 #, tmp218, tmp218, stmp_sum_7.153
	vaddss	%xmm2, %xmm0, %xmm0	 # stmp_sum_7.153, stmp_sum_7.153, stmp_sum_7.153
	vaddss	%xmm1, %xmm0, %xmm0	 # stmp_sum_7.153, stmp_sum_7.153, <retval>
	cmpq	%r9, %r8	 # niters_vector_mult_vf.148, niters.146
	je	.L24	 #,
.L13:
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	leaq	1(%rax), %r9	 #, i
 # sum_vecto_test.cpp:37: 		sum += array[i];
	vaddss	(%rcx,%rax,4), %xmm0, %xmm0	 # *_3, <retval>, <retval>
 # sum_vecto_test.cpp:37: 		sum += array[i];
	leaq	0(,%rax,4), %r8	 #, _2
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	cmpq	%r9, %rdx	 # i, N
	jbe	.L24	 #,
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	leaq	2(%rax), %r9	 #, i
 # sum_vecto_test.cpp:37: 		sum += array[i];
	vaddss	4(%rcx,%r8), %xmm0, %xmm0	 # *_67, <retval>, <retval>
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	cmpq	%r9, %rdx	 # i, N
	jbe	.L24	 #,
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	leaq	3(%rax), %r9	 #, i
 # sum_vecto_test.cpp:37: 		sum += array[i];
	vaddss	8(%rcx,%r8), %xmm0, %xmm0	 # *_113, <retval>, <retval>
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	cmpq	%r9, %rdx	 # i, N
	jbe	.L24	 #,
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	leaq	4(%rax), %r9	 #, i
 # sum_vecto_test.cpp:37: 		sum += array[i];
	vaddss	12(%rcx,%r8), %xmm0, %xmm0	 # *_120, <retval>, <retval>
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	cmpq	%r9, %rdx	 # i, N
	jbe	.L24	 #,
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	leaq	5(%rax), %r9	 #, i
 # sum_vecto_test.cpp:37: 		sum += array[i];
	vaddss	16(%rcx,%r8), %xmm0, %xmm0	 # *_127, <retval>, <retval>
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	cmpq	%r9, %rdx	 # i, N
	jbe	.L24	 #,
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	addq	$6, %rax	 #, i
 # sum_vecto_test.cpp:37: 		sum += array[i];
	vaddss	20(%rcx,%r8), %xmm0, %xmm0	 # *_134, <retval>, <retval>
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	cmpq	%rax, %rdx	 # i, N
	jbe	.L24	 #,
 # sum_vecto_test.cpp:37: 		sum += array[i];
	vaddss	24(%rcx,%r8), %xmm0, %xmm0	 # *_72, <retval>, <retval>
	vzeroupper
 # sum_vecto_test.cpp:40: }
	ret	
	.p2align 4,,10
	.p2align 3
.L24:
	vzeroupper
	ret	
	.p2align 4,,10
	.p2align 3
.L15:
 # sum_vecto_test.cpp:35: 	float sum = 0;
	vxorps	%xmm0, %xmm0, %xmm0	 # <retval>
 # sum_vecto_test.cpp:40: }
	ret	
.L16:
 # sum_vecto_test.cpp:36: 	for(size_t i=0;i<N;++i){
	xorl	%eax, %eax	 # tmp.149
 # sum_vecto_test.cpp:35: 	float sum = 0;
	vxorps	%xmm0, %xmm0, %xmm0	 # <retval>
	jmp	.L9	 #
	.seh_endproc
	.p2align 4
	.globl	_Z6next_ff
	.def	_Z6next_ff;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z6next_ff
_Z6next_ff:
.LFB13478:
	.seh_endprologue
 # sum_vecto_test.cpp:44: 	float y = x + 3.14;
	vcvtss2sd	%xmm0, %xmm0, %xmm0	 # tmp92, tmp88
 # sum_vecto_test.cpp:44: 	float y = x + 3.14;
	vaddsd	.LC1(%rip), %xmm0, %xmm0	 #, tmp88, tmp89
 # sum_vecto_test.cpp:44: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp89, y
 # sum_vecto_test.cpp:45: 	return y - floor(y);
	vrndscaless	$9, %xmm0, %xmm0, %xmm1	 #, y, _6
 # sum_vecto_test.cpp:45: 	return y - floor(y);
	vsubss	%xmm1, %xmm0, %xmm0	 # _6, y, tmp91
 # sum_vecto_test.cpp:46: }
	ret	
	.seh_endproc
	.p2align 4
	.globl	_Z3rndv
	.def	_Z3rndv;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z3rndv
_Z3rndv:
.LFB13479:
	.seh_endprologue
	vxorps	%xmm0, %xmm0, %xmm0	 # tmp92
 # sum_vecto_test.cpp:44: 	float y = x + 3.14;
	vcvtss2sd	_ZL9rand_init(%rip), %xmm0, %xmm0	 # rand_init, tmp92, tmp93
 # sum_vecto_test.cpp:44: 	float y = x + 3.14;
	vaddsd	.LC1(%rip), %xmm0, %xmm0	 #, tmp88, tmp89
 # sum_vecto_test.cpp:44: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp89, y
 # sum_vecto_test.cpp:45: 	return y - floor(y);
	vrndscaless	$9, %xmm0, %xmm0, %xmm1	 #, y, _7
	vsubss	%xmm1, %xmm0, %xmm0	 # _7, y, <retval>
 # sum_vecto_test.cpp:50: 	rand_init = next_f(rand_init);
	vmovss	%xmm0, _ZL9rand_init(%rip)	 # <retval>, rand_init
 # sum_vecto_test.cpp:52: }
	ret	
	.seh_endproc
	.def	__main;	.scl	2;	.type	32;	.endef
	.section .rdata,"dr"
.LC4:
	.ascii "default time = \0"
.LC5:
	.ascii "ms\12\0"
.LC6:
	.ascii "eval simple time = \0"
.LC7:
	.ascii "result = \0"
.LC8:
	.ascii "eval parallel time = \0"
	.section	.text.startup,"x"
	.p2align 4
	.globl	main
	.def	main;	.scl	2;	.type	32;	.endef
	.seh_proc	main
main:
.LFB13480:
	pushq	%r15	 #
	.seh_pushreg	%r15
	pushq	%r14	 #
	.seh_pushreg	%r14
	pushq	%r13	 #
	.seh_pushreg	%r13
	pushq	%r12	 #
	.seh_pushreg	%r12
	pushq	%rsi	 #
	.seh_pushreg	%rsi
	pushq	%rbx	 #
	.seh_pushreg	%rbx
	subq	$104, %rsp	 #,
	.seh_stackalloc	104
	vmovaps	%xmm6, 32(%rsp)	 #,
	.seh_savexmm	%xmm6, 32
	vmovaps	%xmm7, 48(%rsp)	 #,
	.seh_savexmm	%xmm7, 48
	vmovaps	%xmm8, 64(%rsp)	 #,
	.seh_savexmm	%xmm8, 64
	vmovaps	%xmm9, 80(%rsp)	 #,
	.seh_savexmm	%xmm9, 80
	.seh_endprologue
 # sum_vecto_test.cpp:59: 	float * array = (float * ) _aligned_malloc(N*sizeof(float),16*sizeof(float));
	movl	$4000000000, %esi	 #, tmp130
 # sum_vecto_test.cpp:56: int main(int argc,char ** argv){
	call	__main	 #
 # sum_vecto_test.cpp:59: 	float * array = (float * ) _aligned_malloc(N*sizeof(float),16*sizeof(float));
	movl	$64, %edx	 #,
	movq	%rsi, %rcx	 # tmp130,
	call	*__imp__aligned_malloc(%rip)	 #
	movq	%rax, %r12	 # tmp180, _20
 # sum_vecto_test.cpp:61: 	auto t0 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
	vmovss	_ZL9rand_init(%rip), %xmm0	 # rand_init, rand_init_lsm.171
	vmovsd	.LC1(%rip), %xmm2	 #, tmp179
	movq	%rax, %rbx	 # tmp181, t0
	leaq	(%r12,%rsi), %rdx	 #, _16
	movq	%r12, %rax	 # _20, ivtmp.179
	.p2align 4,,10
	.p2align 3
.L28:
 # sum_vecto_test.cpp:44: 	float y = x + 3.14;
	vcvtss2sd	%xmm0, %xmm0, %xmm0	 # rand_init_lsm.171, tmp133
 # sum_vecto_test.cpp:44: 	float y = x + 3.14;
	vaddsd	%xmm2, %xmm0, %xmm0	 # tmp179, tmp133, tmp134
 # sum_vecto_test.cpp:63: 		for(size_t i=0;i<N;++i){
	addq	$4, %rax	 #, ivtmp.179
 # sum_vecto_test.cpp:44: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp134, y
 # sum_vecto_test.cpp:45: 	return y - floor(y);
	vrndscaless	$9, %xmm0, %xmm0, %xmm1	 #, y, _38
	vsubss	%xmm1, %xmm0, %xmm0	 # _38, y, rand_init_lsm.171
 # sum_vecto_test.cpp:64: 			array[i] = rnd();
	vmovss	%xmm0, -4(%rax)	 # rand_init_lsm.171, MEM[(float *)_47]
 # sum_vecto_test.cpp:63: 		for(size_t i=0;i<N;++i){
	cmpq	%rax, %rdx	 # ivtmp.179, _16
	jne	.L28	 #,
	vxorps	%xmm6, %xmm6, %xmm6	 # tmp196
	vmovss	%xmm0, _ZL9rand_init(%rip)	 # rand_init_lsm.171, rand_init
 # sum_vecto_test.cpp:67:     auto t1 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:653: 	return __cd(__cd(__lhs).count() - __cd(__rhs).count());
	subq	%rbx, %rax	 # t0, tmp136
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vcvtsi2ssq	%rax, %xmm6, %xmm0	 # tmp136, tmp196, tmp197
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vmovss	.LC2(%rip), %xmm7	 #, tmp139
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmovss	.LC3(%rip), %xmm8	 #, tmp141
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movq	.refptr._ZSt4cout(%rip), %r13	 #, tmp143
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vdivss	%xmm7, %xmm0, %xmm0	 # tmp139, tmp137, tmp138
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$15, %r8d	 #,
	movq	%r13, %rcx	 # tmp143,
	leaq	.LC4(%rip), %rdx	 #, tmp142
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmulss	%xmm8, %xmm0, %xmm0	 # tmp141, tmp138, tmp140
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:226: 	    return _ToDur(static_cast<__to_rep>(
	vcvttss2siq	%xmm0, %r14	 # tmp140, _43
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r13, %rcx	 # tmp143,
	movq	%r14, %rdx	 # _43,
	call	_ZNSo9_M_insertIxEERSoT_	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	leaq	.LC5(%rip), %r14	 #, tmp145
	movl	$3, %r8d	 #,
	movq	%r14, %rdx	 # tmp145,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%rax, %rcx	 # tmp183, _49
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # sum_vecto_test.cpp:77: 	t0 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # sum_vecto_test.cpp:79: 		result = sum(array,N);
	movq	%r12, %rcx	 # _20,
	movl	$1000000000, %edx	 #,
 # sum_vecto_test.cpp:77: 	t0 = Time::now();
	movq	%rax, %rbx	 # tmp184, t0
 # sum_vecto_test.cpp:79: 		result = sum(array,N);
	call	_Z3sumPKfy	 #
	vmovaps	%xmm0, %xmm9	 # tmp185, result
 # sum_vecto_test.cpp:81: 	t1 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:653: 	return __cd(__cd(__lhs).count() - __cd(__rhs).count());
	subq	%rbx, %rax	 # t0, tmp146
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vcvtsi2ssq	%rax, %xmm6, %xmm2	 # tmp146, tmp196, tmp198
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$19, %r8d	 #,
	movq	%r13, %rcx	 # tmp143,
	leaq	.LC6(%rip), %rdx	 #, tmp152
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vdivss	%xmm7, %xmm2, %xmm2	 # tmp139, tmp147, tmp148
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmulss	%xmm8, %xmm2, %xmm2	 # tmp141, tmp148, tmp150
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:226: 	    return _ToDur(static_cast<__to_rep>(
	vcvttss2siq	%xmm2, %r15	 # tmp150, _51
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r13, %rcx	 # tmp143,
	movq	%r15, %rdx	 # _51,
	call	_ZNSo9_M_insertIxEERSoT_	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	leaq	.LC7(%rip), %r15	 #, tmp156
	movq	%r14, %rdx	 # tmp145,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%rax, %rcx	 # tmp187, _55
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$3, %r8d	 #,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
	movl	$9, %r8d	 #,
	movq	%r15, %rdx	 # tmp156,
	movq	%r13, %rcx	 # tmp143,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # sum_vecto_test.cpp:85: 	std::cout << "result = " << result/N <<std::endl << std::endl;
	vdivss	%xmm7, %xmm9, %xmm1	 # tmp139, result, tmp158
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r13, %rcx	 # tmp143,
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # tmp158,
	call	_ZNSo9_M_insertIdEERSoT_	 #
	movq	%rax, %rcx	 # tmp188, _58
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
	movq	%rax, %rcx	 # tmp189, _59
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
 # sum_vecto_test.cpp:89: 	t0 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # sum_vecto_test.cpp:91: 		result = sum_simd(array,N);
	movq	%r12, %rcx	 # _20,
	movl	$1000000000, %edx	 #,
 # sum_vecto_test.cpp:89: 	t0 = Time::now();
	movq	%rax, %rbx	 # tmp190, t0
 # sum_vecto_test.cpp:91: 		result = sum_simd(array,N);
	call	_Z8sum_simdPKfy	 #
	vmovaps	%xmm0, %xmm9	 # tmp191, result
 # sum_vecto_test.cpp:93: 	t1 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:653: 	return __cd(__cd(__lhs).count() - __cd(__rhs).count());
	subq	%rbx, %rax	 # t0, tmp162
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vcvtsi2ssq	%rax, %xmm6, %xmm6	 # tmp162, tmp196, tmp199
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$21, %r8d	 #,
	movq	%r13, %rcx	 # tmp143,
	leaq	.LC8(%rip), %rdx	 #, tmp168
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vdivss	%xmm7, %xmm6, %xmm2	 # tmp139, tmp163, tmp164
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmulss	%xmm8, %xmm2, %xmm2	 # tmp141, tmp164, tmp166
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:226: 	    return _ToDur(static_cast<__to_rep>(
	vcvttss2siq	%xmm2, %r12	 # tmp166, _62
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r13, %rcx	 # tmp143,
	movq	%r12, %rdx	 # _62,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp193, _66
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movq	%r14, %rdx	 # tmp145,
	movl	$3, %r8d	 #,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
	movl	$9, %r8d	 #,
	movq	%r15, %rdx	 # tmp156,
	movq	%r13, %rcx	 # tmp143,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # sum_vecto_test.cpp:97: 	std::cout << "result = " << result/N <<std::endl << std::endl;
	vdivss	%xmm7, %xmm9, %xmm1	 # tmp139, result, tmp174
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r13, %rcx	 # tmp143,
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # tmp174,
	call	_ZNSo9_M_insertIdEERSoT_	 #
	movq	%rax, %rcx	 # tmp194, _68
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
	movq	%rax, %rcx	 # tmp195, _69
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
	nop	
 # sum_vecto_test.cpp:102: }
	vmovaps	32(%rsp), %xmm6	 #,
	vmovaps	48(%rsp), %xmm7	 #,
	vmovaps	64(%rsp), %xmm8	 #,
	vmovaps	80(%rsp), %xmm9	 #,
	xorl	%eax, %eax	 #
	addq	$104, %rsp	 #,
	popq	%rbx	 #
	popq	%rsi	 #
	popq	%r12	 #
	popq	%r13	 #
	popq	%r14	 #
	popq	%r15	 #
	ret	
	.seh_endproc
	.p2align 4
	.def	_GLOBAL__sub_I__Z8sum_simdPKfy;	.scl	3;	.type	32;	.endef
	.seh_proc	_GLOBAL__sub_I__Z8sum_simdPKfy
_GLOBAL__sub_I__Z8sum_simdPKfy:
.LFB14361:
	subq	$40, %rsp	 #,
	.seh_stackalloc	40
	.seh_endprologue
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/iostream:74:   static ios_base::Init __ioinit;
	leaq	_ZStL8__ioinit(%rip), %rcx	 #, tmp82
	call	_ZNSt8ios_base4InitC1Ev	 #
	leaq	__tcf_1(%rip), %rcx	 #, tmp83
 # sum_vecto_test.cpp:102: }
	addq	$40, %rsp	 #,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/iostream:74:   static ios_base::Init __ioinit;
	jmp	atexit	 #
	.seh_endproc
	.section	.ctors,"w"
	.align 8
	.quad	_GLOBAL__sub_I__Z8sum_simdPKfy
	.data
	.align 4
_ZL9rand_init:
	.long	1094954516
.lcomm _ZStL8__ioinit,1,1
	.section .rdata,"dr"
	.align 8
.LC1:
	.long	1374389535
	.long	1074339512
	.align 4
.LC2:
	.long	1315859240
	.align 4
.LC3:
	.long	1148846080
	.ident	"GCC: (x86_64-posix-seh-rev3, Built by MinGW-W64 project) 11.2.0"
	.def	_ZNSt8ios_base4InitD1Ev;	.scl	2;	.type	32;	.endef
	.def	_ZNSt6chrono3_V212system_clock3nowEv;	.scl	2;	.type	32;	.endef
	.def	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x;	.scl	2;	.type	32;	.endef
	.def	_ZNSo9_M_insertIxEERSoT_;	.scl	2;	.type	32;	.endef
	.def	_ZNSo9_M_insertIdEERSoT_;	.scl	2;	.type	32;	.endef
	.def	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_;	.scl	2;	.type	32;	.endef
	.def	_ZNSt8ios_base4InitC1Ev;	.scl	2;	.type	32;	.endef
	.def	atexit;	.scl	2;	.type	32;	.endef
	.section	.rdata$.refptr._ZSt4cout, "dr"
	.globl	.refptr._ZSt4cout
	.linkonce	discard
.refptr._ZSt4cout:
	.quad	_ZSt4cout
