	.file	"poly_test.cpp"
 # GNU C++23 (x86_64-posix-seh-rev3, Built by MinGW-W64 project) version 11.2.0 (x86_64-w64-mingw32)
 #	compiled by GNU C version 11.2.0, GMP version 6.2.1, MPFR version 4.1.0, MPC version 1.2.1, isl version isl-0.24-GMP

 # GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
 # options passed: -march=alderlake -mmmx -mpopcnt -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx2 -mno-sse4a -mno-fma4 -mno-xop -mfma -mno-avx512f -mbmi -mbmi2 -maes -mpclmul -mno-avx512vl -mno-avx512bw -mno-avx512dq -mno-avx512cd -mno-avx512er -mno-avx512pf -mno-avx512vbmi -mno-avx512ifma -mno-avx5124vnniw -mno-avx5124fmaps -mno-avx512vpopcntdq -mno-avx512vbmi2 -mgfni -mvpclmulqdq -mno-avx512vnni -mno-avx512bitalg -mno-avx512bf16 -mno-avx512vp2intersect -mno-3dnow -madx -mabm -mno-cldemote -mclflushopt -mclwb -mno-clzero -mcx16 -mno-enqcmd -mf16c -mfsgsbase -mfxsr -mno-hle -msahf -mno-lwp -mlzcnt -mmovbe -mmovdir64b -mmovdiri -mno-mwaitx -mno-pconfig -mno-pku -mno-prefetchwt1 -mprfchw -mptwrite -mrdpid -mrdrnd -mrdseed -mno-rtm -mserialize -mno-sgx -msha -mshstk -mno-tbm -mno-tsxldtrk -mvaes -mwaitpkg -mno-wbnoinvd -mxsave -mxsavec -mxsaveopt -mxsaves -mno-amx-tile -mno-amx-int8 -mno-amx-bf16 -mno-uintr -mhreset -mno-kl -mno-widekl -mavxvnni --param=l1-cache-size=48 --param=l1-cache-line-size=64 --param=l2-cache-size=24576 -mtune=alderlake -mavx -O3 -std=c++23
	.text
	.p2align 4
	.def	__tcf_1;	.scl	3;	.type	32;	.endef
	.seh_proc	__tcf_1
__tcf_1:
.LFB14512:
	.seh_endprologue
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/iostream:74:   static ios_base::Init __ioinit;
	leaq	_ZStL8__ioinit(%rip), %rcx	 #, tmp82
	jmp	_ZNSt8ios_base4InitD1Ev	 #
	.seh_endproc
	.p2align 4
	.def	_Z15get_rand_numberPN4evdm15PolynomEstrin12IfEEf.constprop.0;	.scl	3;	.type	32;	.endef
	.seh_proc	_Z15get_rand_numberPN4evdm15PolynomEstrin12IfEEf.constprop.0
_Z15get_rand_numberPN4evdm15PolynomEstrin12IfEEf.constprop.0:
.LFB14516:
	.seh_endprologue
 # poly_test.cpp:147: }
	ret	
	.seh_endproc
	.p2align 4
	.globl	_Z6next_ff
	.def	_Z6next_ff;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z6next_ff
_Z6next_ff:
.LFB13504:
	.seh_endprologue
 # poly_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm0, %xmm0, %xmm0	 # tmp92, tmp88
 # poly_test.cpp:109: 	float y = x + 3.14;
	vaddsd	.LC0(%rip), %xmm0, %xmm0	 #, tmp88, tmp89
 # poly_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp89, y
 # poly_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm1	 #, y, _6
 # poly_test.cpp:110: 	return y - floor(y);
	vsubss	%xmm1, %xmm0, %xmm0	 # _6, y, tmp91
 # poly_test.cpp:111: }
	ret	
	.seh_endproc
	.p2align 4
	.globl	_Z3rndv
	.def	_Z3rndv;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z3rndv
_Z3rndv:
.LFB13505:
	.seh_endprologue
	vxorps	%xmm0, %xmm0, %xmm0	 # tmp92
 # poly_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	_ZL9rand_init(%rip), %xmm0, %xmm0	 # rand_init, tmp92, tmp93
 # poly_test.cpp:109: 	float y = x + 3.14;
	vaddsd	.LC0(%rip), %xmm0, %xmm0	 #, tmp88, tmp89
 # poly_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp89, y
 # poly_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm1	 #, y, _7
	vsubss	%xmm1, %xmm0, %xmm0	 # _7, y, <retval>
 # poly_test.cpp:115: 	rand_init = next_f(rand_init);
	vmovss	%xmm0, _ZL9rand_init(%rip)	 # <retval>, rand_init
 # poly_test.cpp:117: }
	ret	
	.seh_endproc
	.p2align 4
	.globl	_Z15get_rand_numberPN4evdm15PolynomEstrin12IfEEf
	.def	_Z15get_rand_numberPN4evdm15PolynomEstrin12IfEEf;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z15get_rand_numberPN4evdm15PolynomEstrin12IfEEf
_Z15get_rand_numberPN4evdm15PolynomEstrin12IfEEf:
.LFB13508:
	.seh_endprologue
 # poly_test.cpp:145: float  __attribute__ ((noinline)) get_rand_number(poly_t * Poly,float x){
	vmovaps	%xmm1, %xmm0	 # tmp86, x
 # poly_test.cpp:147: }
	ret	
	.seh_endproc
	.p2align 4
	.globl	_Z11eval_hornerPN4evdm15PolynomEstrin12IfEEf
	.def	_Z11eval_hornerPN4evdm15PolynomEstrin12IfEEf;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z11eval_hornerPN4evdm15PolynomEstrin12IfEEf
_Z11eval_hornerPN4evdm15PolynomEstrin12IfEEf:
.LFB13509:
	.seh_endprologue
 # ../src/utils/polynom.hpp:20:             T sum = coeffs[_size-2] + coeffs[_size-1]*x;
	vmovss	44(%rcx), %xmm2	 # MEM[(const float[12] &)Poly_2(D)][11], sum
 # ../src/utils/polynom.hpp:23:                 sum += coeffs[i];
	vmovaps	%xmm1, %xmm0	 # x, x
 # ../src/utils/polynom.hpp:20:             T sum = coeffs[_size-2] + coeffs[_size-1]*x;
	vfmadd213ss	40(%rcx), %xmm1, %xmm2	 # MEM[(const float[12] &)Poly_2(D)][10], x, sum
 # ../src/utils/polynom.hpp:23:                 sum += coeffs[i];
	vfmadd213ss	36(%rcx), %xmm1, %xmm2	 # MEM[(const float[12] &)Poly_2(D)][9], x, sum
	vfmadd213ss	32(%rcx), %xmm1, %xmm2	 # MEM[(const float[12] &)Poly_2(D)][8], x, sum
	vfmadd213ss	28(%rcx), %xmm1, %xmm2	 # MEM[(const float[12] &)Poly_2(D)][7], x, sum
	vfmadd213ss	24(%rcx), %xmm1, %xmm2	 # MEM[(const float[12] &)Poly_2(D)][6], x, sum
	vfmadd213ss	20(%rcx), %xmm1, %xmm2	 # MEM[(const float[12] &)Poly_2(D)][5], x, sum
	vfmadd213ss	16(%rcx), %xmm1, %xmm2	 # MEM[(const float[12] &)Poly_2(D)][4], x, sum
	vfmadd213ss	12(%rcx), %xmm1, %xmm2	 # MEM[(const float[12] &)Poly_2(D)][3], x, sum
	vfmadd213ss	8(%rcx), %xmm1, %xmm2	 # MEM[(const float[12] &)Poly_2(D)][2], x, sum
	vfmadd213ss	4(%rcx), %xmm1, %xmm2	 # MEM[(const float[12] &)Poly_2(D)][1], x, sum
	vfmadd213ss	(%rcx), %xmm2, %xmm0	 # MEM[(const float[12] &)Poly_2(D)][0], sum, x
 # poly_test.cpp:150: }
	ret	
	.seh_endproc
	.section	.text$_ZNK4evdm15PolynomEstrin12IfEclEf,"x"
	.linkonce discard
	.align 2
	.p2align 4
	.globl	_ZNK4evdm15PolynomEstrin12IfEclEf
	.def	_ZNK4evdm15PolynomEstrin12IfEclEf;	.scl	2;	.type	32;	.endef
	.seh_proc	_ZNK4evdm15PolynomEstrin12IfEclEf
_ZNK4evdm15PolynomEstrin12IfEclEf:
.LFB13512:
	.seh_endprologue
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd_builtin.h:1651:       { return __x._M_data + __y._M_data; }
	vmovaps	(%rcx), %xmm4	 # MEM <__int128 unsigned> [(char * {ref-all})_49], tmp116
 # ../src/utils/polynom.hpp:145:             auto y = x*x;
	vmulss	%xmm1, %xmm1, %xmm0	 # x, x, _3
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd.h:1670:       [](auto... __xx) { return __vector_type_t<_Tp, _Np>{__xx...}; },
	vshufps	$0, %xmm1, %xmm1, %xmm3	 # x, tmp106
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd_builtin.h:1651:       { return __x._M_data + __y._M_data; }
	vfmadd132ps	16(%rcx), %xmm4, %xmm3	 # MEM <__int128 unsigned> [(char * {ref-all})_42], tmp116, _52
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd.h:1670:       [](auto... __xx) { return __vector_type_t<_Tp, _Np>{__xx...}; },
	vshufps	$0, %xmm0, %xmm0, %xmm2	 # _3, tmp107
 # ../src/utils/polynom.hpp:148:             y*=x;
	vmulss	%xmm1, %xmm0, %xmm1	 # x, _3, _5
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd_builtin.h:1651:       { return __x._M_data + __y._M_data; }
	vfmadd132ps	32(%rcx), %xmm3, %xmm2	 # MEM <__int128 unsigned> [(char * {ref-all})_18], _52, _25
 # ../src/utils/polynom.hpp:20:             T sum = coeffs[_size-2] + coeffs[_size-1]*x;
	vunpckhps	%xmm2, %xmm2, %xmm3	 # _25, _25, tmp109
	vshufps	$255, %xmm2, %xmm2, %xmm0	 #, _25, _25, tmp108
	vfmadd132ss	%xmm1, %xmm3, %xmm0	 # _5, tmp109, sum
 # ../src/utils/polynom.hpp:23:                 sum += coeffs[i];
	vshufps	$85, %xmm2, %xmm2, %xmm3	 #, _25, _25, tmp110
	vfmadd132ss	%xmm1, %xmm3, %xmm0	 # _5, tmp110, sum
	vfmadd132ss	%xmm1, %xmm2, %xmm0	 # _5, tmp111, <retval>
 # ../src/utils/polynom.hpp:150:         }
	ret	
	.seh_endproc
	.text
	.p2align 4
	.globl	_Z12eval_defaultPN4evdm15PolynomEstrin12IfEEf
	.def	_Z12eval_defaultPN4evdm15PolynomEstrin12IfEEf;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z12eval_defaultPN4evdm15PolynomEstrin12IfEEf
_Z12eval_defaultPN4evdm15PolynomEstrin12IfEEf:
.LFB13511:
	.seh_endprologue
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd_builtin.h:1651:       { return __x._M_data + __y._M_data; }
	vmovaps	(%rcx), %xmm4	 # MEM <__int128 unsigned> [(char * {ref-all})_35], tmp116
 # ../src/utils/polynom.hpp:145:             auto y = x*x;
	vmulss	%xmm1, %xmm1, %xmm3	 # x, x, _8
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd.h:1670:       [](auto... __xx) { return __vector_type_t<_Tp, _Np>{__xx...}; },
	vshufps	$0, %xmm1, %xmm1, %xmm2	 # x, tmp106
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd_builtin.h:1651:       { return __x._M_data + __y._M_data; }
	vfmadd132ps	16(%rcx), %xmm4, %xmm2	 # MEM <__int128 unsigned> [(char * {ref-all})_28], tmp116, _38
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd.h:1670:       [](auto... __xx) { return __vector_type_t<_Tp, _Np>{__xx...}; },
	vshufps	$0, %xmm3, %xmm3, %xmm0	 # _8, tmp107
 # ../src/utils/polynom.hpp:148:             y*=x;
	vmulss	%xmm3, %xmm1, %xmm1	 # _8, x, _20
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd_builtin.h:1651:       { return __x._M_data + __y._M_data; }
	vfmadd132ps	32(%rcx), %xmm2, %xmm0	 # MEM <__int128 unsigned> [(char * {ref-all})_10], _38, _19
 # ../src/utils/polynom.hpp:20:             T sum = coeffs[_size-2] + coeffs[_size-1]*x;
	vunpckhps	%xmm0, %xmm0, %xmm3	 # _19, _19, tmp109
	vshufps	$255, %xmm0, %xmm0, %xmm2	 #, _19, _19, tmp108
	vfmadd132ss	%xmm1, %xmm3, %xmm2	 # _20, tmp109, sum
 # ../src/utils/polynom.hpp:23:                 sum += coeffs[i];
	vshufps	$85, %xmm0, %xmm0, %xmm3	 #, _19, _19, tmp110
	vfmadd132ss	%xmm1, %xmm3, %xmm2	 # _20, tmp110, sum
	vfmadd231ss	%xmm2, %xmm1, %xmm0	 # sum, _20, <retval>
 # poly_test.cpp:153: }
	ret	
	.seh_endproc
	.def	__main;	.scl	2;	.type	32;	.endef
	.section .rdata,"dr"
.LC3:
	.ascii "Coeffs: \0"
.LC4:
	.ascii "Polynom: \0"
.LC5:
	.ascii "P(0) = \0"
.LC6:
	.ascii " vs \0"
.LC7:
	.ascii "P(1) = \0"
.LC9:
	.ascii "P(2) = \0"
.LC13:
	.ascii "default time = \0"
.LC14:
	.ascii "ms\12\0"
.LC15:
	.ascii "sum = \0"
.LC17:
	.ascii "Horner Scheme = \0"
.LC18:
	.ascii "delta time = \0"
.LC19:
	.ascii "Estring Scheme = \0"
.LC20:
	.ascii ", \0"
.LC21:
	.ascii "\0"
	.section	.text.startup,"x"
	.p2align 4
	.globl	main
	.def	main;	.scl	2;	.type	32;	.endef
	.seh_proc	main
main:
.LFB13530:
	pushq	%r15	 #
	.seh_pushreg	%r15
	pushq	%r14	 #
	.seh_pushreg	%r14
	pushq	%r13	 #
	.seh_pushreg	%r13
	pushq	%r12	 #
	.seh_pushreg	%r12
	pushq	%rdi	 #
	.seh_pushreg	%rdi
	pushq	%rsi	 #
	.seh_pushreg	%rsi
	pushq	%rbx	 #
	.seh_pushreg	%rbx
	subq	$336, %rsp	 #,
	.seh_stackalloc	336
	vmovaps	%xmm6, 176(%rsp)	 #,
	.seh_savexmm	%xmm6, 176
	vmovaps	%xmm7, 192(%rsp)	 #,
	.seh_savexmm	%xmm7, 192
	vmovaps	%xmm8, 208(%rsp)	 #,
	.seh_savexmm	%xmm8, 208
	vmovaps	%xmm9, 224(%rsp)	 #,
	.seh_savexmm	%xmm9, 224
	vmovaps	%xmm10, 240(%rsp)	 #,
	.seh_savexmm	%xmm10, 240
	vmovaps	%xmm11, 256(%rsp)	 #,
	.seh_savexmm	%xmm11, 256
	vmovaps	%xmm12, 272(%rsp)	 #,
	.seh_savexmm	%xmm12, 272
	vmovaps	%xmm13, 288(%rsp)	 #,
	.seh_savexmm	%xmm13, 288
	vmovaps	%xmm14, 304(%rsp)	 #,
	.seh_savexmm	%xmm14, 304
	vmovaps	%xmm15, 320(%rsp)	 #,
	.seh_savexmm	%xmm15, 320
	.seh_endprologue
	vxorps	%xmm7, %xmm7, %xmm7	 # tmp491
	leaq	.LC20(%rip), %rbx	 #, tmp448
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	leaq	.LC21(%rip), %rdi	 #, tmp456
 # poly_test.cpp:157: int main(int argc,char ** argv){
	call	__main	 #
 # D:/Soft/Qt/Tools/mingw1120_64/x86_64-w64-mingw32/include/time.h:246: static __inline time_t __CRTDECL time(time_t *_Time) { return _time64(_Time); }
	xorl	%ecx, %ecx	 #
	call	*__imp__time64(%rip)	 #
 # poly_test.cpp:168: 	srand(time(0));
	movl	%eax, %ecx	 # tmp457, _57
	call	srand	 #
 # poly_test.cpp:170: 	rand_init = rand()/(RAND_MAX+0.0f);
	call	rand	 #
 # poly_test.cpp:170: 	rand_init = rand()/(RAND_MAX+0.0f);
	vcvtsi2ssl	%eax, %xmm7, %xmm0	 # tmp458, tmp491, tmp492
 # poly_test.cpp:173: 	cout << "Coeffs: ";
	movq	.refptr._ZSt4cout(%rip), %r12	 #, tmp449
 # poly_test.cpp:157: int main(int argc,char ** argv){
	leaq	127(%rsp), %rsi	 #, tmp267
 # poly_test.cpp:173: 	cout << "Coeffs: ";
	leaq	.LC3(%rip), %rdx	 #, tmp278
	movq	%r12, %rcx	 # tmp449,
 # poly_test.cpp:170: 	rand_init = rand()/(RAND_MAX+0.0f);
	vdivss	.LC2(%rip), %xmm0, %xmm0	 #, tmp275, tmp276
 # poly_test.cpp:157: int main(int argc,char ** argv){
	andq	$-32, %rsi	 #, tmp269
 # poly_test.cpp:170: 	rand_init = rand()/(RAND_MAX+0.0f);
	vmovss	%xmm0, _ZL9rand_init(%rip)	 # tmp276, rand_init
 # poly_test.cpp:173: 	cout << "Coeffs: ";
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc	 #
	vmovsd	.LC0(%rip), %xmm6	 #, tmp450
 # poly_test.cpp:174: 	for(size_t j = 0;j<16;++j){
	xorl	%eax, %eax	 # j
	jmp	.L20	 #
	.p2align 4,,10
	.p2align 3
.L11:
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movq	%r12, %rcx	 # tmp449,
	movl	$2, %r8d	 #,
	movq	%rbx, %rdx	 # tmp448,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp449,
	vcvtss2sd	%xmm8, %xmm8, %xmm1	 # _44,
	call	_ZNSo9_M_insertIdEERSoT_	 #
 # poly_test.cpp:174: 	for(size_t j = 0;j<16;++j){
	cmpq	$16, %r13	 #, _366
	je	.L19	 #,
.L24:
 # poly_test.cpp:174: 	for(size_t j = 0;j<16;++j){
	movq	%r13, %rax	 # _366, j
.L20:
 # poly_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	_ZL9rand_init(%rip), %xmm7, %xmm8	 # rand_init, tmp491, tmp493
 # poly_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm8, %xmm8	 # tmp450, tmp280, tmp281
 # poly_test.cpp:174: 	for(size_t j = 0;j<16;++j){
	leaq	1(%rax), %r13	 #, _366
 # poly_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm8, %xmm8, %xmm8	 # tmp281, y
 # poly_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm8, %xmm8, %xmm0	 #, y, _46
	vsubss	%xmm0, %xmm8, %xmm8	 # _46, y, _44
 # poly_test.cpp:115: 	rand_init = next_f(rand_init);
	vmovss	%xmm8, _ZL9rand_init(%rip)	 # _44, rand_init
 # poly_test.cpp:175: 		coeffs[j] = rnd();
	vmovss	%xmm8, (%rsi,%rax,4)	 # _44, MEM[(float *)&coeffs + j_304 * 4]
 # poly_test.cpp:176: 		cout << (j ? ", " : "") << coeffs[j];
	testq	%rax, %rax	 # j
	jne	.L11	 #,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movq	%r12, %rcx	 # tmp449,
	xorl	%r8d, %r8d	 #
	movq	%rdi, %rdx	 # tmp456,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp449,
	vcvtss2sd	%xmm8, %xmm8, %xmm1	 # _44,
	call	_ZNSo9_M_insertIdEERSoT_	 #
	jmp	.L24	 #
.L19:
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	movq	%r12, %rcx	 # tmp449,
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
 # ../src/utils/polynom.hpp:160:             for(size_t i=0;i<N;++i){
	xorl	%edi, %edi	 # i
	.p2align 4,,10
	.p2align 3
.L13:
 # ../src/utils/polynom.hpp:137:             auto dv = std::div(i,3);
	movl	$3, %edx	 #,
	movl	%edi, %ecx	 # i, i
 # ../src/utils/polynom.hpp:161:                 (*this)[i] = cfs[i];
	vmovss	(%rsi,%rdi,4), %xmm8	 # MEM[(float *)&coeffs + i_320 * 4], _242
 # ../src/utils/polynom.hpp:137:             auto dv = std::div(i,3);
	call	div	 #
 # ../src/utils/polynom.hpp:138:             return coeffs[ (dv.rem)*4 +dv.quot];
	movq	%rax, %rdx	 # tmp459, tmp285
	shrq	$32, %rdx	 #, tmp285
 # ../src/utils/polynom.hpp:138:             return coeffs[ (dv.rem)*4 +dv.quot];
	leal	(%rax,%rdx,4), %eax	 #, _247
 # ../src/utils/polynom.hpp:161:                 (*this)[i] = cfs[i];
	cltq
 # ../src/utils/polynom.hpp:160:             for(size_t i=0;i<N;++i){
	incq	%rdi	 # i
 # ../src/utils/polynom.hpp:161:                 (*this)[i] = cfs[i];
	vmovss	%xmm8, 48(%rsp,%rax,4)	 # _242, MEM <const struct PolynomEstrin12> [(float &)&P].coeffs[_247]
 # ../src/utils/polynom.hpp:160:             for(size_t i=0;i<N;++i){
	cmpq	$12, %rdi	 #, i
	jne	.L13	 #,
 # poly_test.cpp:179: 	cout << "Polynom: ";
	leaq	.LC4(%rip), %rdx	 #, tmp288
	movq	%r12, %rcx	 # tmp449,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc	 #
 # poly_test.cpp:180: 	for(size_t i=0;i<P.size();++i){
	xorl	%eax, %eax	 # i
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	leaq	.LC21(%rip), %r13	 #, tmp455
 # ../src/utils/polynom.hpp:137:             auto dv = std::div(i,3);
	movl	%eax, %r14d	 # i, _364
 # poly_test.cpp:180: 	for(size_t i=0;i<P.size();++i){
	leaq	1(%rax), %rdi	 #, _365
 # poly_test.cpp:181: 		cout << (i ? ", " : "") << P[i];
	testq	%rax, %rax	 # i
	je	.L33	 #,
.L14:
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$2, %r8d	 #,
	movq	%rbx, %rdx	 # tmp448,
	movq	%r12, %rcx	 # tmp449,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # ../src/utils/polynom.hpp:137:             auto dv = std::div(i,3);
	movl	$3, %edx	 #,
	movl	%r14d, %ecx	 # _364,
	call	div	 #
 # ../src/utils/polynom.hpp:138:             return coeffs[ (dv.rem)*4 +dv.quot];
	movq	%rax, %rdx	 # tmp489, tmp428
	shrq	$32, %rdx	 #, tmp428
 # ../src/utils/polynom.hpp:138:             return coeffs[ (dv.rem)*4 +dv.quot];
	leal	(%rax,%rdx,4), %eax	 #, _85
 # poly_test.cpp:181: 		cout << (i ? ", " : "") << P[i];
	cltq
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp449,
	vcvtss2sd	48(%rsp,%rax,4), %xmm7, %xmm1	 # MEM <struct poly_t> [(float &)&P].coeffs[_85], tmp491, tmp497
	call	_ZNSo9_M_insertIdEERSoT_	 #
 # poly_test.cpp:180: 	for(size_t i=0;i<P.size();++i){
	cmpq	$12, %rdi	 #, _365
	je	.L21	 #,
.L23:
 # poly_test.cpp:180: 	for(size_t i=0;i<P.size();++i){
	movq	%rdi, %rax	 # _365, i
 # ../src/utils/polynom.hpp:137:             auto dv = std::div(i,3);
	movl	%eax, %r14d	 # i, _364
 # poly_test.cpp:180: 	for(size_t i=0;i<P.size();++i){
	leaq	1(%rax), %rdi	 #, _365
 # poly_test.cpp:181: 		cout << (i ? ", " : "") << P[i];
	testq	%rax, %rax	 # i
	jne	.L14	 #,
.L33:
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	xorl	%r8d, %r8d	 #
	movq	%r13, %rdx	 # tmp455,
	movq	%r12, %rcx	 # tmp449,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # ../src/utils/polynom.hpp:137:             auto dv = std::div(i,3);
	movl	%r14d, %ecx	 # _364,
	movl	$3, %edx	 #,
	call	div	 #
 # ../src/utils/polynom.hpp:138:             return coeffs[ (dv.rem)*4 +dv.quot];
	movq	%rax, %rdx	 # tmp490, tmp435
	shrq	$32, %rdx	 #, tmp435
 # ../src/utils/polynom.hpp:138:             return coeffs[ (dv.rem)*4 +dv.quot];
	leal	(%rax,%rdx,4), %eax	 #, _94
 # poly_test.cpp:181: 		cout << (i ? ", " : "") << P[i];
	cltq
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp449,
	vcvtss2sd	48(%rsp,%rax,4), %xmm7, %xmm1	 # MEM <struct poly_t> [(float &)&P].coeffs[_94], tmp491, tmp498
	call	_ZNSo9_M_insertIdEERSoT_	 #
	jmp	.L23	 #
	.p2align 4,,10
	.p2align 3
.L21:
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	movq	%r12, %rcx	 # tmp449,
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
 # poly_test.cpp:183: 	cout << "P(0) = " << P(0) << " vs " << coeff_evaulator1<Deg>::eval(coeffs,0.0f) << endl;
	leaq	48(%rsp), %rbx	 #, tmp446
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$7, %r8d	 #,
	leaq	.LC5(%rip), %rdx	 #, tmp291
	movq	%r12, %rcx	 # tmp449,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_test.cpp:183: 	cout << "P(0) = " << P(0) << " vs " << coeff_evaulator1<Deg>::eval(coeffs,0.0f) << endl;
	movq	%rbx, %rcx	 # tmp446,
	vxorps	%xmm1, %xmm1, %xmm1	 #
	call	_ZNK4evdm15PolynomEstrin12IfEclEf	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp449,
	vcvtss2sd	%xmm0, %xmm0, %xmm1	 # tmp460,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	leaq	.LC6(%rip), %r14	 #, tmp296
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	call	_ZNSo9_M_insertIdEERSoT_	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movq	%rax, %rcx	 # _89,
	movq	%r14, %rdx	 # tmp296,
	movl	$4, %r8d	 #,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%rax, %r13	 # tmp461, _89
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_test.cpp:24: 			sum += data[i];
	vxorps	%xmm11, %xmm11, %xmm11	 # sum
	vaddss	44(%rsi), %xmm11, %xmm8	 # MEM[(const float *)&coeffs + 44B], sum, sum
 # poly_test.cpp:24: 			sum += data[i];
	movl	40(%rsi), %edi	 # MEM[(const float *)&coeffs + 40B], _434
	vmovss	36(%rsi), %xmm5	 # MEM[(const float *)&coeffs + 36B], _445
 # poly_test.cpp:24: 			sum += data[i];
	vmovd	%edi, %xmm1	 # _434, sum
	vfmadd231ss	%xmm11, %xmm8, %xmm1	 # sum, sum, sum
 # poly_test.cpp:24: 			sum += data[i];
	vmovss	32(%rsi), %xmm4	 # MEM[(const float *)&coeffs + 32B], _456
 # poly_test.cpp:24: 			sum += data[i];
	vmovss	%xmm5, 44(%rsp)	 # _445, %sfp
 # poly_test.cpp:24: 			sum += data[i];
	movl	24(%rsi), %r15d	 # MEM[(const float *)&coeffs + 24B], _478
	vmovss	20(%rsi), %xmm15	 # MEM[(const float *)&coeffs + 20B], _489
 # poly_test.cpp:24: 			sum += data[i];
	vfmadd132ss	%xmm11, %xmm5, %xmm1	 # sum, _445, sum
 # poly_test.cpp:24: 			sum += data[i];
	vmovss	28(%rsi), %xmm5	 # MEM[(const float *)&coeffs + 28B], _467
 # poly_test.cpp:24: 			sum += data[i];
	vmovd	%r15d, %xmm3	 # _478, _478
 # poly_test.cpp:24: 			sum += data[i];
	vmovss	16(%rsi), %xmm14	 # MEM[(const float *)&coeffs + 16B], _500
	vmovss	12(%rsi), %xmm13	 # MEM[(const float *)&coeffs + 12B], _511
 # poly_test.cpp:24: 			sum += data[i];
	vfmadd132ss	%xmm11, %xmm4, %xmm1	 # sum, _456, sum
 # poly_test.cpp:24: 			sum += data[i];
	vmovss	8(%rsi), %xmm12	 # MEM[(const float *)&coeffs + 8B], _522
	vmovss	4(%rsi), %xmm10	 # MEM[(const float *)&coeffs + 4B], _533
	vmovss	(%rsi), %xmm9	 # MEM[(const float *)&coeffs], _98
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r13, %rcx	 # _89,
 # poly_test.cpp:24: 			sum += data[i];
	vfmadd132ss	%xmm11, %xmm5, %xmm1	 # sum, _467, sum
	vmovss	%xmm4, 36(%rsp)	 # _456, %sfp
	vmovss	%xmm5, 40(%rsp)	 # _467, %sfp
	vfmadd132ss	%xmm11, %xmm3, %xmm1	 # sum, _478, sum
	vfmadd132ss	%xmm11, %xmm15, %xmm1	 # sum, _489, sum
	vfmadd132ss	%xmm11, %xmm14, %xmm1	 # sum, _500, sum
	vfmadd132ss	%xmm11, %xmm13, %xmm1	 # sum, _511, sum
	vfmadd132ss	%xmm11, %xmm12, %xmm1	 # sum, _522, sum
	vfmadd132ss	%xmm11, %xmm10, %xmm1	 # sum, _533, sum
	vfmadd132ss	%xmm11, %xmm9, %xmm1	 # sum, _98, sum
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # sum,
	call	_ZNSo9_M_insertIdEERSoT_	 #
	movq	%rax, %rcx	 # tmp462, _91
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$7, %r8d	 #,
	leaq	.LC7(%rip), %rdx	 #, tmp311
	movq	%r12, %rcx	 # tmp449,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_test.cpp:184: 	cout << "P(1) = " << P(1) << " vs " << coeff_evaulator1<Deg>::eval(coeffs,1.0f) << endl;
	vmovss	.LC8(%rip), %xmm1	 #,
	movq	%rbx, %rcx	 # tmp446,
	call	_ZNK4evdm15PolynomEstrin12IfEclEf	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp449,
	vcvtss2sd	%xmm0, %xmm0, %xmm1	 # tmp463,
	call	_ZNSo9_M_insertIdEERSoT_	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movq	%rax, %rcx	 # _103,
	movl	$4, %r8d	 #,
	movq	%r14, %rdx	 # tmp296,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%rax, %r13	 # tmp464, _103
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_test.cpp:24: 			sum += data[i];
	vmovups	16(%rsi), %ymm1	 # MEM <const vector(8) float> [(const float *)&coeffs + 16B], MEM <const vector(8) float> [(const float *)&coeffs + 16B]
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r13, %rcx	 # _103,
	vextractf128	$0x1, %ymm1, %xmm0	 # MEM <const vector(8) float> [(const float *)&coeffs + 16B], tmp318
	vshufps	$255, %xmm0, %xmm0, %xmm2	 #, tmp318, tmp318, stmp_sum_114.231
	vaddss	%xmm11, %xmm2, %xmm2	 # sum, stmp_sum_114.231, stmp_sum_114.231
	vunpckhps	%xmm0, %xmm0, %xmm3	 # tmp318, tmp318, stmp_sum_114.231
	vshufps	$85, %xmm1, %xmm1, %xmm4	 #, tmp323, tmp323, stmp_sum_114.231
	vaddss	%xmm2, %xmm3, %xmm2	 # stmp_sum_114.231, stmp_sum_114.231, stmp_sum_114.231
	vshufps	$85, %xmm0, %xmm0, %xmm3	 #, tmp318, tmp318, stmp_sum_114.231
	vaddss	%xmm2, %xmm3, %xmm2	 # stmp_sum_114.231, stmp_sum_114.231, stmp_sum_114.231
	vshufps	$255, %xmm1, %xmm1, %xmm3	 #, tmp323, tmp323, stmp_sum_114.231
	vaddss	%xmm2, %xmm0, %xmm0	 # stmp_sum_114.231, stmp_sum_114.231, stmp_sum_114.231
	vaddss	%xmm0, %xmm3, %xmm3	 # stmp_sum_114.231, stmp_sum_114.231, stmp_sum_114.231
	vunpckhps	%xmm1, %xmm1, %xmm0	 # tmp323, tmp323, stmp_sum_114.231
	vaddss	%xmm3, %xmm0, %xmm3	 # stmp_sum_114.231, stmp_sum_114.231, stmp_sum_114.231
 # poly_test.cpp:24: 			sum += data[i];
	vmovaps	(%rsi), %xmm0	 # MEM <const vector(4) float> [(const float *)&coeffs], vect__326.234
 # poly_test.cpp:24: 			sum += data[i];
	vaddss	%xmm3, %xmm4, %xmm3	 # stmp_sum_114.231, stmp_sum_114.231, stmp_sum_114.231
	vshufps	$255, %xmm0, %xmm0, %xmm2	 #, vect__326.234, vect__326.234, stmp_sum_325.236
	vaddss	%xmm1, %xmm3, %xmm1	 # stmp_sum_114.231, stmp_sum_114.231, sum
	vaddss	%xmm2, %xmm1, %xmm1	 # stmp_sum_325.236, sum, stmp_sum_325.236
	vunpckhps	%xmm0, %xmm0, %xmm2	 # vect__326.234, vect__326.234, stmp_sum_325.236
	vaddss	%xmm1, %xmm2, %xmm2	 # stmp_sum_325.236, stmp_sum_325.236, stmp_sum_325.236
	vshufps	$85, %xmm0, %xmm0, %xmm1	 #, vect__326.234, vect__326.234, stmp_sum_325.236
	vaddss	%xmm2, %xmm1, %xmm1	 # stmp_sum_325.236, stmp_sum_325.236, stmp_sum_325.236
	vaddss	%xmm0, %xmm1, %xmm1	 # stmp_sum_325.236, stmp_sum_325.236, sum
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # sum,
	vzeroupper
	call	_ZNSo9_M_insertIdEERSoT_	 #
	movq	%rax, %rcx	 # tmp465, _106
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$7, %r8d	 #,
	leaq	.LC9(%rip), %rdx	 #, tmp334
	movq	%r12, %rcx	 # tmp449,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_test.cpp:185: 	cout << "P(2) = " << P(2) << " vs " << coeff_evaulator1<Deg>::eval(coeffs,2.0f) << endl;
	vmovss	.LC10(%rip), %xmm1	 #,
	movq	%rbx, %rcx	 # tmp446,
	call	_ZNK4evdm15PolynomEstrin12IfEclEf	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp449,
	vcvtss2sd	%xmm0, %xmm0, %xmm1	 # tmp466,
	call	_ZNSo9_M_insertIdEERSoT_	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movq	%rax, %rcx	 # _117,
	movl	$4, %r8d	 #,
	movq	%r14, %rdx	 # tmp296,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%rax, %r13	 # tmp467, _117
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_test.cpp:24: 			sum += data[i];
	vmovss	.LC10(%rip), %xmm0	 #, tmp340
	vmovd	%edi, %xmm1	 # _434, _434
	vfmadd231ss	%xmm0, %xmm8, %xmm1	 # tmp340, sum, _434
	vmovss	44(%rsp), %xmm5	 # %sfp, _445
	vmovd	%r15d, %xmm4	 # _478, _478
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r13, %rcx	 # _117,
 # poly_test.cpp:24: 			sum += data[i];
	vfmadd132ss	%xmm0, %xmm5, %xmm1	 # tmp340, _445, sum
	vfmadd213ss	36(%rsp), %xmm0, %xmm1	 # %sfp, tmp340, sum
	vfmadd213ss	40(%rsp), %xmm0, %xmm1	 # %sfp, tmp340, sum
	vfmadd132ss	%xmm0, %xmm4, %xmm1	 # tmp340, _478, sum
	vfmadd132ss	%xmm0, %xmm15, %xmm1	 # tmp340, _489, sum
	vfmadd132ss	%xmm0, %xmm14, %xmm1	 # tmp340, _500, sum
	vfmadd132ss	%xmm0, %xmm13, %xmm1	 # tmp340, _511, sum
	vfmadd132ss	%xmm0, %xmm12, %xmm1	 # tmp340, _522, sum
	vfmadd132ss	%xmm0, %xmm10, %xmm1	 # tmp340, _533, sum
	vfmadd132ss	%xmm0, %xmm9, %xmm1	 # tmp340, _98, sum
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # sum,
	call	_ZNSo9_M_insertIdEERSoT_	 #
	movq	%rax, %rcx	 # tmp468, _119
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
 # poly_test.cpp:187: 	auto t0 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
	vmovss	_ZL9rand_init(%rip), %xmm1	 # rand_init, rand_init_lsm.223
	movq	%rax, %rsi	 # tmp469, t0
 # poly_test.cpp:167: 	float sum = 0;
	vmovaps	%xmm11, %xmm12	 # sum, sum
 # poly_test.cpp:187: 	auto t0 = Time::now();
	movl	$100000000, %eax	 #, ivtmp_347
	.p2align 4,,10
	.p2align 3
.L16:
 # poly_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # rand_init_lsm.223, tmp352
 # poly_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm1, %xmm1	 # tmp450, tmp352, tmp353
 # poly_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm1, %xmm1, %xmm1	 # tmp353, y
 # poly_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm1, %xmm1, %xmm0	 #, y, _133
	vsubss	%xmm0, %xmm1, %xmm1	 # _133, y, rand_init_lsm.223
 # poly_test.cpp:194: 		sum += get_rand_number(&P,rnd());
	vmovaps	%xmm1, %xmm0	 # rand_init_lsm.223,
	call	_Z15get_rand_numberPN4evdm15PolynomEstrin12IfEEf.constprop.0	 #
 # poly_test.cpp:194: 		sum += get_rand_number(&P,rnd());
	vaddss	%xmm0, %xmm12, %xmm12	 # tmp470, sum, sum
 # poly_test.cpp:188: 	for(size_t i=0;i<Ntest;++i){
	decq	%rax	 # ivtmp_347
	jne	.L16	 #,
	vmovss	%xmm1, _ZL9rand_init(%rip)	 # rand_init_lsm.223, rand_init
 # poly_test.cpp:196:     auto t1 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:653: 	return __cd(__cd(__lhs).count() - __cd(__rhs).count());
	subq	%rsi, %rax	 # t0, tmp355
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vcvtsi2ssq	%rax, %xmm7, %xmm0	 # tmp355, tmp491, tmp494
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vmovss	.LC11(%rip), %xmm10	 #, tmp444
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmovss	.LC12(%rip), %xmm9	 #, tmp445
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$15, %r8d	 #,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vdivss	%xmm10, %xmm0, %xmm0	 # tmp444, tmp356, tmp357
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	leaq	.LC13(%rip), %rdx	 #, tmp361
	movq	%r12, %rcx	 # tmp449,
 # poly_test.cpp:199:     std::cout << "default time = " << d.count() << "ms\n";
	leaq	.LC14(%rip), %r13	 #, tmp452
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	leaq	.LC15(%rip), %r14	 #, tmp453
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmulss	%xmm9, %xmm0, %xmm0	 # tmp445, tmp357, tmp359
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:226: 	    return _ToDur(static_cast<__to_rep>(
	vcvttss2siq	%xmm0, %rsi	 # tmp359, _136
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r12, %rcx	 # tmp449,
	movq	%rsi, %rdx	 # _136,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp472, _143
 # poly_test.cpp:199:     std::cout << "default time = " << d.count() << "ms\n";
	movq	%r13, %rdx	 # tmp452,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$6, %r8d	 #,
	movq	%r14, %rdx	 # tmp453,
	movq	%r12, %rcx	 # tmp449,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_test.cpp:200: 	std::cout << "sum = " << sum/Ntest <<std::endl << std::endl;
	vmovss	.LC16(%rip), %xmm8	 #, tmp447
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp449,
 # poly_test.cpp:200: 	std::cout << "sum = " << sum/Ntest <<std::endl << std::endl;
	vdivss	%xmm8, %xmm12, %xmm1	 # tmp447, sum, tmp367
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # tmp367,
	call	_ZNSo9_M_insertIdEERSoT_	 #
	movq	%rax, %rcx	 # tmp473, _145
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
	movq	%rax, %rcx	 # tmp474, _146
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
 # poly_test.cpp:205: 	t0 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # poly_test.cpp:115: 	rand_init = next_f(rand_init);
	vmovss	_ZL9rand_init(%rip), %xmm1	 # rand_init, _152
 # poly_test.cpp:205: 	t0 = Time::now();
	movq	%rax, %rdi	 # tmp475, t0
 # poly_test.cpp:204: 	sum = 0;
	vmovaps	%xmm11, %xmm12	 # sum, sum
 # poly_test.cpp:115: 	rand_init = next_f(rand_init);
	movl	$100000000, %eax	 #, ivtmp_349
	.p2align 4,,10
	.p2align 3
.L17:
 # poly_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # _152, tmp371
 # poly_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm1, %xmm1	 # tmp450, tmp371, tmp372
 # poly_test.cpp:212: 		sum += eval_horner(&P,rnd());
	movq	%rbx, %rcx	 # tmp446,
 # poly_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm1, %xmm1, %xmm1	 # tmp372, y
 # poly_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm1, %xmm1, %xmm0	 #, y, _151
	vsubss	%xmm0, %xmm1, %xmm1	 # _151, y, _152
 # poly_test.cpp:115: 	rand_init = next_f(rand_init);
	vmovss	%xmm1, _ZL9rand_init(%rip)	 # _152, rand_init
 # poly_test.cpp:212: 		sum += eval_horner(&P,rnd());
	call	_Z11eval_hornerPN4evdm15PolynomEstrin12IfEEf	 #
 # poly_test.cpp:212: 		sum += eval_horner(&P,rnd());
	vaddss	%xmm0, %xmm12, %xmm12	 # tmp476, sum, sum
 # poly_test.cpp:206: 	for(size_t i=0;i<Ntest;++i){
	decq	%rax	 # ivtmp_349
	jne	.L17	 #,
 # poly_test.cpp:214: 	t1 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:653: 	return __cd(__cd(__lhs).count() - __cd(__rhs).count());
	subq	%rdi, %rax	 # t0, tmp375
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vcvtsi2ssq	%rax, %xmm7, %xmm0	 # tmp375, tmp491, tmp495
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$16, %r8d	 #,
	leaq	.LC17(%rip), %rdx	 #, tmp381
	movq	%r12, %rcx	 # tmp449,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vdivss	%xmm10, %xmm0, %xmm0	 # tmp444, tmp376, tmp377
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	leaq	.LC18(%rip), %r15	 #, tmp454
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmulss	%xmm9, %xmm0, %xmm0	 # tmp445, tmp377, tmp379
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:226: 	    return _ToDur(static_cast<__to_rep>(
	vcvttss2siq	%xmm0, %rdi	 # tmp379, _154
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r12, %rcx	 # tmp449,
	movq	%rdi, %rdx	 # _154,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp478, _158
 # poly_test.cpp:217: 	std::cout << "Horner Scheme = " << d.count() << "ms\n";
	movq	%r13, %rdx	 # tmp452,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc	 #
 # poly_test.cpp:218: 	std::cout << "delta time = " << d.count() - def_time  << "ms\n";
	subq	%rsi, %rdi	 # _136, _154
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$13, %r8d	 #,
	movq	%r15, %rdx	 # tmp454,
	movq	%r12, %rcx	 # tmp449,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_test.cpp:218: 	std::cout << "delta time = " << d.count() - def_time  << "ms\n";
	movq	%rdi, %rdx	 # _154, tmp387
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r12, %rcx	 # tmp449,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp479, _159
 # poly_test.cpp:218: 	std::cout << "delta time = " << d.count() - def_time  << "ms\n";
	movq	%r13, %rdx	 # tmp452,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$6, %r8d	 #,
	movq	%r14, %rdx	 # tmp453,
	movq	%r12, %rcx	 # tmp449,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_test.cpp:219: 	std::cout << "sum = " << sum/Ntest << std::endl << std::endl;
	vdivss	%xmm8, %xmm12, %xmm1	 # tmp447, sum, tmp392
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp449,
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # tmp392,
	call	_ZNSo9_M_insertIdEERSoT_	 #
	movq	%rax, %rcx	 # tmp480, _161
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
	movq	%rax, %rcx	 # tmp481, _162
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
 # poly_test.cpp:222: 	t0 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # poly_test.cpp:115: 	rand_init = next_f(rand_init);
	vmovss	_ZL9rand_init(%rip), %xmm1	 # rand_init, _168
 # poly_test.cpp:222: 	t0 = Time::now();
	movq	%rax, %rdi	 # tmp482, t0
 # poly_test.cpp:115: 	rand_init = next_f(rand_init);
	movl	$100000000, %eax	 #, ivtmp_351
	.p2align 4,,10
	.p2align 3
.L18:
 # poly_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # _168, tmp396
 # poly_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm1, %xmm1	 # tmp450, tmp396, tmp397
 # poly_test.cpp:229: 		sum += eval_default(&P,rnd());
	movq	%rbx, %rcx	 # tmp446,
 # poly_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm1, %xmm1, %xmm1	 # tmp397, y
 # poly_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm1, %xmm1, %xmm0	 #, y, _167
	vsubss	%xmm0, %xmm1, %xmm1	 # _167, y, _168
 # poly_test.cpp:115: 	rand_init = next_f(rand_init);
	vmovss	%xmm1, _ZL9rand_init(%rip)	 # _168, rand_init
 # poly_test.cpp:229: 		sum += eval_default(&P,rnd());
	vmovss	%xmm1, 36(%rsp)	 # _168, %sfp
	call	_Z12eval_defaultPN4evdm15PolynomEstrin12IfEEf	 #
 # poly_test.cpp:223: 	for(size_t i=0;i<Ntest;++i){
	decq	%rax	 # ivtmp_351
 # poly_test.cpp:229: 		sum += eval_default(&P,rnd());
	vaddss	%xmm0, %xmm11, %xmm11	 # tmp483, sum, sum
 # poly_test.cpp:223: 	for(size_t i=0;i<Ntest;++i){
	vmovss	36(%rsp), %xmm1	 # %sfp, _168
	jne	.L18	 #,
 # poly_test.cpp:231: 	t1 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:653: 	return __cd(__cd(__lhs).count() - __cd(__rhs).count());
	subq	%rdi, %rax	 # t0, tmp400
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vcvtsi2ssq	%rax, %xmm7, %xmm7	 # tmp400, tmp491, tmp496
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$17, %r8d	 #,
	movq	%r12, %rcx	 # tmp449,
	leaq	.LC19(%rip), %rdx	 #, tmp406
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vdivss	%xmm10, %xmm7, %xmm0	 # tmp444, tmp401, tmp402
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmulss	%xmm9, %xmm0, %xmm0	 # tmp445, tmp402, tmp404
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:226: 	    return _ToDur(static_cast<__to_rep>(
	vcvttss2siq	%xmm0, %rbx	 # tmp404, _170
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r12, %rcx	 # tmp449,
	movq	%rbx, %rdx	 # _170,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp485, _174
 # poly_test.cpp:234: 	std::cout << "Estring Scheme = " << d.count() << "ms\n";
	movq	%r13, %rdx	 # tmp452,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$13, %r8d	 #,
	movq	%r15, %rdx	 # tmp454,
	movq	%r12, %rcx	 # tmp449,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_test.cpp:235: 	std::cout << "delta time = " << d.count() - def_time  << "ms\n";
	movq	%rbx, %rdx	 # _170, _170
	subq	%rsi, %rdx	 # _136, _170
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r12, %rcx	 # tmp449,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp486, _175
 # poly_test.cpp:235: 	std::cout << "delta time = " << d.count() - def_time  << "ms\n";
	movq	%r13, %rdx	 # tmp452,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$6, %r8d	 #,
	movq	%r14, %rdx	 # tmp453,
	movq	%r12, %rcx	 # tmp449,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_test.cpp:236: 	std::cout << "sum = " << sum/Ntest << std::endl << std::endl;
	vdivss	%xmm8, %xmm11, %xmm1	 # tmp447, sum, tmp417
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp449,
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # tmp417,
	call	_ZNSo9_M_insertIdEERSoT_	 #
	movq	%rax, %rcx	 # tmp487, _177
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
	movq	%rax, %rcx	 # tmp488, _178
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_	 #
	nop	
 # poly_test.cpp:240: }
	vmovaps	176(%rsp), %xmm6	 #,
	vmovaps	192(%rsp), %xmm7	 #,
	vmovaps	208(%rsp), %xmm8	 #,
	vmovaps	224(%rsp), %xmm9	 #,
	vmovaps	240(%rsp), %xmm10	 #,
	vmovaps	256(%rsp), %xmm11	 #,
	vmovaps	272(%rsp), %xmm12	 #,
	vmovaps	288(%rsp), %xmm13	 #,
	vmovaps	304(%rsp), %xmm14	 #,
	vmovaps	320(%rsp), %xmm15	 #,
	xorl	%eax, %eax	 #
	addq	$336, %rsp	 #,
	popq	%rbx	 #
	popq	%rsi	 #
	popq	%rdi	 #
	popq	%r12	 #
	popq	%r13	 #
	popq	%r14	 #
	popq	%r15	 #
	ret	
	.seh_endproc
	.p2align 4
	.def	_GLOBAL__sub_I__Z6next_ff;	.scl	3;	.type	32;	.endef
	.seh_proc	_GLOBAL__sub_I__Z6next_ff
_GLOBAL__sub_I__Z6next_ff:
.LFB14513:
	subq	$40, %rsp	 #,
	.seh_stackalloc	40
	.seh_endprologue
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/iostream:74:   static ios_base::Init __ioinit;
	leaq	_ZStL8__ioinit(%rip), %rcx	 #, tmp82
	call	_ZNSt8ios_base4InitC1Ev	 #
	leaq	__tcf_1(%rip), %rcx	 #, tmp83
 # poly_test.cpp:240: }
	addq	$40, %rsp	 #,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/iostream:74:   static ios_base::Init __ioinit;
	jmp	atexit	 #
	.seh_endproc
	.section	.ctors,"w"
	.align 8
	.quad	_GLOBAL__sub_I__Z6next_ff
	.data
	.align 4
_ZL9rand_init:
	.long	1094954516
.lcomm _ZStL8__ioinit,1,1
	.section .rdata,"dr"
	.align 8
.LC0:
	.long	1374389535
	.long	1074339512
	.align 4
.LC2:
	.long	1191181824
	.align 4
.LC8:
	.long	1065353216
	.align 4
.LC10:
	.long	1073741824
	.align 4
.LC11:
	.long	1315859240
	.align 4
.LC12:
	.long	1148846080
	.align 4
.LC16:
	.long	1287568416
	.ident	"GCC: (x86_64-posix-seh-rev3, Built by MinGW-W64 project) 11.2.0"
	.def	_ZNSt8ios_base4InitD1Ev;	.scl	2;	.type	32;	.endef
	.def	srand;	.scl	2;	.type	32;	.endef
	.def	rand;	.scl	2;	.type	32;	.endef
	.def	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc;	.scl	2;	.type	32;	.endef
	.def	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x;	.scl	2;	.type	32;	.endef
	.def	_ZNSo9_M_insertIdEERSoT_;	.scl	2;	.type	32;	.endef
	.def	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_;	.scl	2;	.type	32;	.endef
	.def	div;	.scl	2;	.type	32;	.endef
	.def	_ZNSt6chrono3_V212system_clock3nowEv;	.scl	2;	.type	32;	.endef
	.def	_ZNSo9_M_insertIxEERSoT_;	.scl	2;	.type	32;	.endef
	.def	_ZNSt8ios_base4InitC1Ev;	.scl	2;	.type	32;	.endef
	.def	atexit;	.scl	2;	.type	32;	.endef
	.section	.rdata$.refptr._ZSt4cout, "dr"
	.globl	.refptr._ZSt4cout
	.linkonce	discard
.refptr._ZSt4cout:
	.quad	_ZSt4cout
