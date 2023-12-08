	.file	"poly_eval_test.cpp"
 # GNU C++23 (x86_64-posix-seh-rev3, Built by MinGW-W64 project) version 11.2.0 (x86_64-w64-mingw32)
 #	compiled by GNU C version 11.2.0, GMP version 6.2.1, MPFR version 4.1.0, MPC version 1.2.1, isl version isl-0.24-GMP

 # GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
 # options passed: -march=alderlake -mmmx -mpopcnt -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx2 -mno-sse4a -mno-fma4 -mno-xop -mfma -mno-avx512f -mbmi -mbmi2 -maes -mpclmul -mno-avx512vl -mno-avx512bw -mno-avx512dq -mno-avx512cd -mno-avx512er -mno-avx512pf -mno-avx512vbmi -mno-avx512ifma -mno-avx5124vnniw -mno-avx5124fmaps -mno-avx512vpopcntdq -mno-avx512vbmi2 -mgfni -mvpclmulqdq -mno-avx512vnni -mno-avx512bitalg -mno-avx512bf16 -mno-avx512vp2intersect -mno-3dnow -madx -mabm -mno-cldemote -mclflushopt -mclwb -mno-clzero -mcx16 -mno-enqcmd -mf16c -mfsgsbase -mfxsr -mno-hle -msahf -mno-lwp -mlzcnt -mmovbe -mmovdir64b -mmovdiri -mno-mwaitx -mno-pconfig -mno-pku -mno-prefetchwt1 -mprfchw -mptwrite -mrdpid -mrdrnd -mrdseed -mno-rtm -mserialize -mno-sgx -msha -mshstk -mno-tbm -mno-tsxldtrk -mvaes -mwaitpkg -mno-wbnoinvd -mxsave -mxsavec -mxsaveopt -mxsaves -mno-amx-tile -mno-amx-int8 -mno-amx-bf16 -mno-uintr -mhreset -mno-kl -mno-widekl -mavxvnni --param=l1-cache-size=48 --param=l1-cache-line-size=64 --param=l2-cache-size=24576 -mtune=alderlake -mavx -O3 -std=c++23
	.text
	.p2align 4
	.def	__tcf_1;	.scl	3;	.type	32;	.endef
	.seh_proc	__tcf_1
__tcf_1:
.LFB14745:
	.seh_endprologue
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/iostream:74:   static ios_base::Init __ioinit;
	leaq	_ZStL8__ioinit(%rip), %rcx	 #, tmp82
	jmp	_ZNSt8ios_base4InitD1Ev	 #
	.seh_endproc
	.p2align 4
	.def	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0;	.scl	3;	.type	32;	.endef
	.seh_proc	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0
_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0:
.LFB14749:
	pushq	%r13	 #
	.seh_pushreg	%r13
	pushq	%r12	 #
	.seh_pushreg	%r12
	subq	$40, %rsp	 #,
	.seh_stackalloc	40
	.seh_endprologue
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:682:     { return flush(__os.put(__os.widen('\n'))); }
	movq	(%rcx), %rax	 # __os_1(D)->_vptr.basic_ostream, __os_1(D)->_vptr.basic_ostream
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:681:     endl(basic_ostream<_CharT, _Traits>& __os)
	movq	%rcx, %r12	 # tmp99, __os
	movq	-24(%rax), %rax	 # MEM[(long long int *)_2 + -24B], MEM[(long long int *)_2 + -24B]
	movq	240(%rcx,%rax), %r13	 # MEM[(const struct __ctype_type * *)_5 + 240B], _14
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/bits/basic_ios.h:49:       if (!__f)
	testq	%r13, %r13	 # _14
	je	.L7	 #,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/bits/locale_facets.h:877: 	if (_M_widen_ok)
	cmpb	$0, 56(%r13)	 #, MEM[(const struct ctype *)_14]._M_widen_ok
	je	.L5	 #,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/bits/locale_facets.h:878: 	  return _M_widen[static_cast<unsigned char>(__c)];
	movsbl	67(%r13), %edx	 # MEM[(const struct ctype *)_14]._M_widen[10],
.L6:
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:682:     { return flush(__os.put(__os.widen('\n'))); }
	movq	%r12, %rcx	 # __os,
	call	_ZNSo3putEc	 #
	movq	%rax, %rcx	 # tmp101, _8
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:682:     { return flush(__os.put(__os.widen('\n'))); }
	addq	$40, %rsp	 #,
	popq	%r12	 #
	popq	%r13	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:704:     { return __os.flush(); }
	jmp	_ZNSo5flushEv	 #
.L5:
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/bits/locale_facets.h:879: 	this->_M_widen_init();
	movq	%r13, %rcx	 # _14,
	call	_ZNKSt5ctypeIcE13_M_widen_initEv	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/bits/locale_facets.h:880: 	return this->do_widen(__c);
	movq	0(%r13), %rax	 # MEM[(const struct ctype *)_14].D.56597._vptr.facet, MEM[(const struct ctype *)_14].D.56597._vptr.facet
	movl	$10, %edx	 #,
	movq	%r13, %rcx	 # _14,
	call	*48(%rax)	 # MEM[(int (*) () *)_24 + 48B]
	movsbl	%al, %edx	 # tmp100,
	jmp	.L6	 #
.L7:
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/bits/basic_ios.h:50: 	__throw_bad_cast();
	call	_ZSt16__throw_bad_castv	 #
	nop	
	.seh_endproc
	.p2align 4
	.def	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0;	.scl	3;	.type	32;	.endef
	.seh_proc	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0:
.LFB14751:
	pushq	%r13	 #
	.seh_pushreg	%r13
	pushq	%r12	 #
	.seh_pushreg	%r12
	subq	$40, %rsp	 #,
	.seh_stackalloc	40
	.seh_endprologue
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:606:     operator<<(basic_ostream<char, _Traits>& __out, const char* __s)
	movq	%rcx, %r13	 # tmp99, __out
	movq	%rdx, %r12	 # tmp100, __s
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:608:       if (!__s)
	testq	%rdx, %rdx	 # __s
	je	.L10	 #,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/bits/char_traits.h:371: 	return __builtin_strlen(__s);
	movq	%rdx, %rcx	 # __s,
	call	strlen	 #
	movq	%rax, %r8	 #, tmp101
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movq	%r12, %rdx	 # __s,
	movq	%r13, %rcx	 # __out,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:614:     }
	addq	$40, %rsp	 #,
	popq	%r12	 #
	popq	%r13	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	jmp	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
.L10:
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:609: 	__out.setstate(ios_base::badbit);
	movq	(%rcx), %rax	 # __out_2(D)->_vptr.basic_ostream, __out_2(D)->_vptr.basic_ostream
	movq	-24(%rax), %rcx	 # MEM[(long long int *)_9 + -24B], __out
	addq	%r13, %rcx	 # __out, __out
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/bits/ios_base.h:170:   { return _Ios_Iostate(static_cast<int>(__a) | static_cast<int>(__b)); }
	movl	32(%rcx), %edx	 # MEM[(const struct basic_ios *)_12].D.60127._M_streambuf_state, tmp94
	orl	$1, %edx	 #, tmp94
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:614:     }
	addq	$40, %rsp	 #,
	popq	%r12	 #
	popq	%r13	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/bits/basic_ios.h:158:       { this->clear(this->rdstate() | __state); }
	jmp	_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate	 #
	.seh_endproc
	.p2align 4
	.globl	_Z6next_ff
	.def	_Z6next_ff;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z6next_ff
_Z6next_ff:
.LFB13479:
	.seh_endprologue
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm0, %xmm0, %xmm0	 # tmp92, tmp88
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	.LC0(%rip), %xmm0, %xmm0	 #, tmp88, tmp89
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp89, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm1	 #, y, _6
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vsubss	%xmm1, %xmm0, %xmm0	 # _6, y, tmp91
 # poly_eval_test.cpp:111: }
	ret	
	.seh_endproc
	.p2align 4
	.globl	_Z3rndv
	.def	_Z3rndv;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z3rndv
_Z3rndv:
.LFB13480:
	.seh_endprologue
	vxorps	%xmm0, %xmm0, %xmm0	 # tmp92
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	_ZL9rand_init(%rip), %xmm0, %xmm0	 # rand_init, tmp92, tmp93
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	.LC0(%rip), %xmm0, %xmm0	 #, tmp88, tmp89
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp89, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm1	 #, y, _7
	vsubss	%xmm1, %xmm0, %xmm0	 # _7, y, <retval>
 # poly_eval_test.cpp:115: 	rand_init = next_f(rand_init);
	vmovss	%xmm0, _ZL9rand_init(%rip)	 # <retval>, rand_init
 # poly_eval_test.cpp:117: }
	ret	
	.seh_endproc
	.section	.text$_Z4evalILy8EEfPff,"x"
	.linkonce discard
	.p2align 4
	.globl	_Z4evalILy8EEfPff
	.def	_Z4evalILy8EEfPff;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z4evalILy8EEfPff
_Z4evalILy8EEfPff:
.LFB13960:
	.seh_endprologue
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd.h:1670:       [](auto... __xx) { return __vector_type_t<_Tp, _Np>{__xx...}; },
	vshufps	$0, %xmm1, %xmm1, %xmm0	 # x, tmp102
 # poly_eval_test.cpp:36:         return eval2deg_poly(std::get<0>(result_tp),std::get<1>(result_tp),x*x);
	vmulss	%xmm1, %xmm1, %xmm1	 # x, x, _48
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd_builtin.h:1651:       { return __x._M_data + __y._M_data; }
	vmovaps	(%rcx), %xmm4	 # MEM <__int128 unsigned> [(char * {ref-all})_9], tmp114
	vfmadd132ps	16(%rcx), %xmm4, %xmm0	 # MEM <__int128 unsigned> [(char * {ref-all})_9 + 16B], tmp114, _53
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd.h:1670:       [](auto... __xx) { return __vector_type_t<_Tp, _Np>{__xx...}; },
	vmovsldup	%xmm1, %xmm2	 # _48, tmp104
 # poly_eval_test.cpp:36:         return eval2deg_poly(std::get<0>(result_tp),std::get<1>(result_tp),x*x);
	vmulss	%xmm1, %xmm1, %xmm1	 # _48, _48, tmp107
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/xmmintrin.h:791:   return (__m128) __builtin_ia32_movhlps ((__v4sf)__A, (__v4sf)__B);
	vmovhlps	%xmm0, %xmm0, %xmm3	 # _53, _53, tmp103
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd_builtin.h:1651:       { return __x._M_data + __y._M_data; }
	vfmadd231ps	%xmm3, %xmm2, %xmm0	 # tmp105, tmp104, _63
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd_scalar.h:212: 			      + __promote_preserving_unsigned(__y));
	vmovshdup	%xmm0, %xmm2	 # _63, tmp108
	vfmadd231ss	%xmm2, %xmm1, %xmm0	 # tmp108, tmp107, <retval>
 # poly_eval_test.cpp:76: }
	ret	
	.seh_endproc
	.section	.text$_Z5eval1ILy8EEfPff,"x"
	.linkonce discard
	.p2align 4
	.globl	_Z5eval1ILy8EEfPff
	.def	_Z5eval1ILy8EEfPff;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z5eval1ILy8EEfPff
_Z5eval1ILy8EEfPff:
.LFB13981:
	.seh_endprologue
 # poly_eval_test.cpp:24: 			sum += data[i];
	vxorps	%xmm2, %xmm2, %xmm2	 # sum
	vfmadd213ss	28(%rcx), %xmm1, %xmm2	 # MEM[(const float *)_3 + 28B], x, sum
	vmovaps	%xmm1, %xmm0	 # x, x
	vfmadd213ss	24(%rcx), %xmm1, %xmm2	 # MEM[(const float *)_3 + 24B], x, sum
	vfmadd213ss	20(%rcx), %xmm1, %xmm2	 # MEM[(const float *)_3 + 20B], x, sum
	vfmadd213ss	16(%rcx), %xmm1, %xmm2	 # MEM[(const float *)_3 + 16B], x, sum
	vfmadd213ss	12(%rcx), %xmm1, %xmm2	 # MEM[(const float *)_3 + 12B], x, sum
	vfmadd213ss	8(%rcx), %xmm1, %xmm2	 # MEM[(const float *)_3 + 8B], x, sum
	vfmadd213ss	4(%rcx), %xmm1, %xmm2	 # MEM[(const float *)_3 + 4B], x, sum
	vfmadd213ss	(%rcx), %xmm2, %xmm0	 # MEM[(const float *)_3], sum, x
 # poly_eval_test.cpp:68: }
	ret	
	.seh_endproc
	.section	.text$_Z5eval2ILy8EEfPff,"x"
	.linkonce discard
	.p2align 4
	.globl	_Z5eval2ILy8EEfPff
	.def	_Z5eval2ILy8EEfPff;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z5eval2ILy8EEfPff
_Z5eval2ILy8EEfPff:
.LFB13982:
	.seh_endprologue
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd_builtin.h:1651:       { return __x._M_data + __y._M_data; }
	vmovaps	(%rcx), %xmm4	 # MEM <__int128 unsigned> [(char * {ref-all})_9], tmp106
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd.h:1670:       [](auto... __xx) { return __vector_type_t<_Tp, _Np>{__xx...}; },
	vshufps	$0, %xmm1, %xmm1, %xmm0	 # x, tmp97
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/experimental/bits/simd_builtin.h:1651:       { return __x._M_data + __y._M_data; }
	vfmadd132ps	16(%rcx), %xmm4, %xmm0	 # MEM <__int128 unsigned> [(char * {ref-all})_9 + 16B], tmp106, _53
 # poly_eval_test.cpp:51: 				sum+=result[i];
	vunpckhps	%xmm0, %xmm0, %xmm3	 # _53, _53, tmp99
	vshufps	$255, %xmm0, %xmm0, %xmm2	 #, _53, _53, tmp98
	vfmadd132ss	%xmm1, %xmm3, %xmm2	 # x, tmp99, sum
	vshufps	$85, %xmm0, %xmm0, %xmm3	 #, _53, _53, tmp100
	vfmadd132ss	%xmm1, %xmm3, %xmm2	 # x, tmp100, sum
	vfmadd231ss	%xmm2, %xmm1, %xmm0	 # sum, x, <retval>
 # poly_eval_test.cpp:106: }
	ret	
	.seh_endproc
	.def	__main;	.scl	2;	.type	32;	.endef
	.section .rdata,"dr"
.LC8:
	.ascii "default time = \0"
.LC9:
	.ascii "ms\12\0"
.LC10:
	.ascii "sum = \0"
.LC12:
	.ascii "eval zero time = \0"
.LC13:
	.ascii "eval parallel time = \0"
.LC14:
	.ascii "delta time = \0"
.LC15:
	.ascii "eval not parallel time = \0"
.LC16:
	.ascii "not simd parallel time = \0"
	.section	.text.startup,"x"
	.p2align 4
	.globl	main
	.def	main;	.scl	2;	.type	32;	.endef
	.seh_proc	main
main:
.LFB13481:
	pushq	%r15	 #
	.seh_pushreg	%r15
	pushq	%r14	 #
	.seh_pushreg	%r14
	pushq	%r13	 #
	.seh_pushreg	%r13
	pushq	%r12	 #
	.seh_pushreg	%r12
	pushq	%rbp	 #
	.seh_pushreg	%rbp
	pushq	%rdi	 #
	.seh_pushreg	%rdi
	pushq	%rsi	 #
	.seh_pushreg	%rsi
	pushq	%rbx	 #
	.seh_pushreg	%rbx
	subq	$760, %rsp	 #,
	.seh_stackalloc	760
	vmovaps	%xmm6, 592(%rsp)	 #,
	.seh_savexmm	%xmm6, 592
	vmovaps	%xmm7, 608(%rsp)	 #,
	.seh_savexmm	%xmm7, 608
	vmovaps	%xmm8, 624(%rsp)	 #,
	.seh_savexmm	%xmm8, 624
	vmovaps	%xmm9, 640(%rsp)	 #,
	.seh_savexmm	%xmm9, 640
	vmovaps	%xmm10, 656(%rsp)	 #,
	.seh_savexmm	%xmm10, 656
	vmovaps	%xmm11, 672(%rsp)	 #,
	.seh_savexmm	%xmm11, 672
	vmovaps	%xmm12, 688(%rsp)	 #,
	.seh_savexmm	%xmm12, 688
	vmovaps	%xmm13, 704(%rsp)	 #,
	.seh_savexmm	%xmm13, 704
	vmovaps	%xmm14, 720(%rsp)	 #,
	.seh_savexmm	%xmm14, 720
	vmovaps	%xmm15, 736(%rsp)	 #,
	.seh_savexmm	%xmm15, 736
	.seh_endprologue
	vxorps	%xmm7, %xmm7, %xmm7	 # tmp548
 # poly_eval_test.cpp:120: int main(int argc,char ** argv){
	leaq	319(%rsp), %rbx	 #, tmp285
	call	__main	 #
 # poly_eval_test.cpp:126: 	alignas(256) float coeffs[Deg+1] = {1,2,3,4,5,6,7,8}; //needs 3 operation op(a,b,x) = a + b*x
	vmovaps	.LC2(%rip), %ymm0	 #, tmp291
 # poly_eval_test.cpp:120: int main(int argc,char ** argv){
	xorb	%bl, %bl	 # tmp287
 # poly_eval_test.cpp:126: 	alignas(256) float coeffs[Deg+1] = {1,2,3,4,5,6,7,8}; //needs 3 operation op(a,b,x) = a + b*x
	vmovaps	%ymm0, 256(%rbx)	 # tmp291, MEM <vector(8) float> [(float *)&coeffs]
 # poly_eval_test.cpp:127: 	alignas(256) std::array<float,8> cf{1,5,3,7,2,6,4,8}; //needs 1 parallel and 1 single operation op(a,b,x)
	vmovaps	.LC3(%rip), %ymm0	 #, tmp292
 # D:/Soft/Qt/Tools/mingw1120_64/x86_64-w64-mingw32/include/time.h:246: static __inline time_t __CRTDECL time(time_t *_Time) { return _time64(_Time); }
	xorl	%ecx, %ecx	 #
 # poly_eval_test.cpp:127: 	alignas(256) std::array<float,8> cf{1,5,3,7,2,6,4,8}; //needs 1 parallel and 1 single operation op(a,b,x)
	vmovaps	%ymm0, (%rbx)	 # tmp292, MEM <vector(8) float> [(float *)&cf]
 # D:/Soft/Qt/Tools/mingw1120_64/x86_64-w64-mingw32/include/time.h:246: static __inline time_t __CRTDECL time(time_t *_Time) { return _time64(_Time); }
	vzeroupper
	call	*__imp__time64(%rip)	 #
 # poly_eval_test.cpp:132: 	srand(time(0));
	movl	%eax, %ecx	 # tmp520, _87
	call	srand	 #
 # poly_eval_test.cpp:134: 	rand_init = rand()/(RAND_MAX+0.0f);
	call	rand	 #
 # poly_eval_test.cpp:134: 	rand_init = rand()/(RAND_MAX+0.0f);
	vcvtsi2ssl	%eax, %xmm7, %xmm3	 # tmp521, tmp548, tmp549
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vmovsd	.LC0(%rip), %xmm6	 #, tmp516
 # poly_eval_test.cpp:134: 	rand_init = rand()/(RAND_MAX+0.0f);
	vdivss	.LC4(%rip), %xmm3, %xmm4	 #, tmp295, _4
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm4, %xmm4, %xmm0	 # _4, tmp297
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp297, tmp298
 # poly_eval_test.cpp:134: 	rand_init = rand()/(RAND_MAX+0.0f);
	vmovd	%xmm4, %esi	 # _4, _4
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp298, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm2	 #, y, _98
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vsubss	%xmm2, %xmm0, %xmm3	 # _98, y, tmp300
	vcvtss2sd	%xmm3, %xmm3, %xmm3	 # tmp300, tmp301
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm3, %xmm3	 # tmp516, tmp301, tmp302
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm3, %xmm3, %xmm3	 # tmp302, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm3, %xmm3, %xmm4	 #, y, _310
	vsubss	%xmm4, %xmm3, %xmm3	 # _310, y, _309
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm3, %xmm3, %xmm0	 # _309, tmp304
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp304, tmp305
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp305, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm5	 #, y, _298
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vsubss	%xmm5, %xmm0, %xmm13	 # _298, y, tmp307
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vmovd	%xmm5, %edx	 # _298, _298
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm13, %xmm13, %xmm13	 # tmp307, tmp308
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm13, %xmm13	 # tmp516, tmp308, tmp309
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm13, %xmm13, %xmm13	 # tmp309, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm13, %xmm13, %xmm14	 #, y, _292
	vsubss	%xmm14, %xmm13, %xmm13	 # _292, y, _291
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm13, %xmm13, %xmm0	 # _291, tmp311
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp311, tmp312
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp312, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm10	 #, y, _280
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vsubss	%xmm10, %xmm0, %xmm0	 # _280, y, tmp314
	vcvtss2sd	%xmm0, %xmm0, %xmm0	 # tmp314, tmp315
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp315, tmp316
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp316, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm8	 #, y, _274
	vsubss	%xmm8, %xmm0, %xmm5	 # _274, y, _273
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm5, %xmm5, %xmm0	 # _273, tmp318
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp318, tmp319
	vmovd	%xmm5, %r9d	 # _273, _273
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp319, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm5	 #, y, _245
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vsubss	%xmm5, %xmm0, %xmm1	 # _245, y, tmp321
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vmovd	%xmm5, %ecx	 # _245, _245
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # tmp321, tmp322
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm1, %xmm1	 # tmp516, tmp322, tmp323
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm1, %xmm1, %xmm1	 # tmp323, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm1, %xmm1, %xmm15	 #, y, _239
	vsubss	%xmm15, %xmm1, %xmm5	 # _239, y, _238
	vunpcklps	%xmm15, %xmm8, %xmm8	 # _239, _274, tmp393
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm5, %xmm5, %xmm0	 # _238, tmp325
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp325, tmp326
	vmovd	%xmm5, %r10d	 # _238, _238
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp326, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm12	 #, y, _225
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vsubss	%xmm12, %xmm0, %xmm0	 # _225, y, tmp328
	vcvtss2sd	%xmm0, %xmm0, %xmm0	 # tmp328, tmp329
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp329, tmp330
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp330, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm5	 #, y, _62
	vsubss	%xmm5, %xmm0, %xmm0	 # _62, y, _26
	vmovd	%xmm0, %r11d	 # _26, _26
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm0, %xmm0, %xmm0	 # _26, tmp332
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp332, tmp333
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp333, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm1	 #, y, _28
	vmovd	%xmm1, %r8d	 # _28, _28
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vsubss	%xmm1, %xmm0, %xmm1	 # _28, y, tmp335
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # tmp335, tmp336
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm1, %xmm1	 # tmp516, tmp336, tmp337
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm1, %xmm1, %xmm1	 # tmp337, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm1, %xmm1, %xmm0	 #, y, _33
	vsubss	%xmm0, %xmm1, %xmm1	 # _33, y, _340
	vmovd	%xmm0, %eax	 # _33, _33
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm1, %xmm1, %xmm0	 # _340, tmp339
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp339, tmp340
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp340, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm11	 #, y, _351
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vsubss	%xmm11, %xmm0, %xmm0	 # _351, y, tmp342
	vcvtss2sd	%xmm0, %xmm0, %xmm0	 # tmp342, tmp343
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp343, tmp344
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp344, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm9	 #, y, _357
	vsubss	%xmm9, %xmm0, %xmm0	 # _357, y, _358
	vmovd	%xmm0, %edi	 # _358, _358
	vunpcklps	%xmm0, %xmm1, %xmm0	 # _358, _340, tmp348
	vmovq	%xmm0, %r14	 # tmp348, tmp348
	vmovd	%r10d, %xmm1	 # _238, _238
	vmovd	%r11d, %xmm0	 # _26, _26
	vunpcklps	%xmm0, %xmm1, %xmm1	 # _26, _238, tmp349
	vmovd	%r9d, %xmm0	 # _273, _273
	vunpcklps	%xmm0, %xmm13, %xmm13	 # _273, _291, tmp351
	vmovd	%esi, %xmm0	 # _4, _4
	vunpcklps	%xmm3, %xmm0, %xmm3	 # _309, _4, tmp352
	vmovlhps	%xmm13, %xmm3, %xmm3	 # tmp351, tmp352, tmp350
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vmovapd	.LC5(%rip), %ymm13	 #, tmp358
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vmovd	%edi, %xmm0	 # _358, _358
	vcvtss2sd	%xmm0, %xmm0, %xmm0	 # _358, tmp353
	vcvtps2pd	%xmm3, %ymm3	 # tmp350, vect__257.472
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp353, tmp354
	vaddpd	%ymm13, %ymm3, %ymm3	 # tmp358, vect__257.472, vect__103.473
	vmovq	%xmm1, %r10	 # tmp349, tmp349
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp354, y
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vmovupd	%ymm3, 32(%rsp)	 # vect__103.473, %sfp
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm1	 #, y, _109
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vmovq	%r14, %xmm3	 # tmp348, tmp348
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vmovd	%xmm1, %r9d	 # _109, _109
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vmovq	%r10, %xmm1	 # tmp349, tmp349
	vmovlhps	%xmm3, %xmm1, %xmm1	 # tmp348, tmp349, tmp360
	vcvtps2pd	%xmm1, %ymm1	 # tmp360, vect__257.472
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddpd	%ymm13, %ymm1, %ymm1	 # tmp358, vect__257.472, vect__103.473
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtpd2psy	32(%rsp), %xmm3	 # %sfp, tmp364
	vcvtpd2psy	%ymm1, %xmm1	 # vect__103.473, tmp365
	vinsertf128	$0x1, %xmm1, %ymm3, %ymm3	 # tmp365, tmp364, vect_y_258.474
	vmovd	%r9d, %xmm1	 # _109, _109
	vunpcklps	%xmm1, %xmm11, %xmm11	 # _109, _351, tmp368
	vmovd	%r8d, %xmm1	 # _28, _28
	vunpcklps	%xmm1, %xmm12, %xmm1	 # _28, _225, tmp369
	vmovlhps	%xmm11, %xmm1, %xmm1	 # tmp368, tmp369, tmp367
	vmovd	%ecx, %xmm11	 # _245, _245
	vunpcklps	%xmm11, %xmm10, %xmm10	 # _245, _280, tmp371
	vmovd	%edx, %xmm11	 # _298, _298
	vunpcklps	%xmm11, %xmm2, %xmm2	 # _298, _98, tmp372
	vmovlhps	%xmm10, %xmm2, %xmm2	 # tmp371, tmp372, tmp370
	vinsertf128	$0x1, %xmm1, %ymm2, %ymm1	 # tmp367, tmp370, tmp366
	vsubps	%ymm1, %ymm3, %ymm1	 # tmp366, vect_y_258.474, vect__259.475
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vmovd	%r9d, %xmm2	 # _109, _109
	vsubss	%xmm2, %xmm0, %xmm0	 # _109, y, tmp373
 # poly_eval_test.cpp:138: 		coeffs[j] = rnd();
	vmovaps	%ymm1, 256(%rbx)	 # vect__259.475, MEM <vector(8) float> [(float *)&coeffs]
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtps2pd	%xmm1, %ymm2	 # vect__259.475, vect__313.478
	vextractf128	$0x1, %ymm1, %xmm1	 # vect__259.475, tmp382
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddpd	%ymm13, %ymm2, %ymm2	 # tmp358, vect__313.478, vect__312.479
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtps2pd	%xmm1, %ymm1	 # tmp382, vect__313.478
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddpd	%ymm13, %ymm1, %ymm1	 # tmp358, vect__313.478, vect__312.479
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm0, %xmm0, %xmm0	 # tmp373, tmp374
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp374, tmp375
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtpd2psy	%ymm2, %xmm2	 # vect__312.479, tmp386
	vcvtpd2psy	%ymm1, %xmm1	 # vect__312.479, tmp387
	vinsertf128	$0x1, %xmm1, %ymm2, %ymm1	 # tmp387, tmp386, vect_y_311.480
	vmovd	%eax, %xmm2	 # _33, _33
	vunpcklps	%xmm2, %xmm5, %xmm2	 # _33, _62, tmp391
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp375, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm0	 #, y, _73
	vunpcklps	%xmm0, %xmm9, %xmm9	 # _73, _357, tmp390
	vunpcklps	%xmm14, %xmm4, %xmm0	 # _292, _310, tmp394
	vmovlhps	%xmm9, %xmm2, %xmm2	 # tmp390, tmp391, tmp389
	vmovlhps	%xmm8, %xmm0, %xmm0	 # tmp393, tmp394, tmp392
	vinsertf128	$0x1, %xmm2, %ymm0, %ymm0	 # tmp389, tmp392, tmp388
	vsubps	%ymm0, %ymm1, %ymm0	 # tmp388, vect_y_311.480, vect__309.481
 # poly_eval_test.cpp:139: 		cf[j] = rnd();
	vmovaps	%ymm0, (%rbx)	 # vect__309.481, MEM <const vector(8) float> [(value_type &)&cf]
	vextractf128	$0x1, %ymm0, %xmm0	 # vect__309.481, tmp396
	vextractps	$3, %xmm0, _ZL9rand_init(%rip)	 #, tmp396, rand_init
 # poly_eval_test.cpp:142: 	auto t0 = Time::now();
	vzeroupper
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # poly_eval_test.cpp:131: 	float sum = 0;
	vxorps	%xmm11, %xmm11, %xmm11	 # sum
	vmovss	_ZL9rand_init(%rip), %xmm0	 # rand_init, rand_init_lsm.464
 # poly_eval_test.cpp:142: 	auto t0 = Time::now();
	movq	%rax, %rsi	 # tmp522, t0
 # poly_eval_test.cpp:131: 	float sum = 0;
	vmovaps	%xmm11, %xmm12	 # sum, sum
 # poly_eval_test.cpp:142: 	auto t0 = Time::now();
	movl	$100000000, %eax	 #, ivtmp_316
	.p2align 4,,10
	.p2align 3
.L17:
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm0, %xmm0, %xmm0	 # rand_init_lsm.464, tmp397
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp397, tmp398
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp398, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm1	 #, y, _117
	vsubss	%xmm1, %xmm0, %xmm0	 # _117, y, rand_init_lsm.464
 # poly_eval_test.cpp:149: 		sum += rnd();
	vaddss	%xmm0, %xmm12, %xmm12	 # rand_init_lsm.464, sum, sum
 # poly_eval_test.cpp:143: 	for(size_t i=0;i<Ntest;++i){
	decq	%rax	 # ivtmp_316
	jne	.L17	 #,
	vmovss	%xmm0, _ZL9rand_init(%rip)	 # rand_init_lsm.464, rand_init
 # poly_eval_test.cpp:151:     auto t1 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:653: 	return __cd(__cd(__lhs).count() - __cd(__rhs).count());
	subq	%rsi, %rax	 # t0, tmp400
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vcvtsi2ssq	%rax, %xmm7, %xmm0	 # tmp400, tmp548, tmp550
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vmovss	.LC6(%rip), %xmm10	 #, tmp518
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmovss	.LC7(%rip), %xmm9	 #, tmp519
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movq	.refptr._ZSt4cout(%rip), %r12	 #, tmp510
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vdivss	%xmm10, %xmm0, %xmm0	 # tmp518, tmp401, tmp402
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$15, %r8d	 #,
	leaq	.LC8(%rip), %rdx	 #, tmp406
	movq	%r12, %rcx	 # tmp510,
 # poly_eval_test.cpp:154:     std::cout << "default time = " << d.count() << "ms\n";
	leaq	.LC9(%rip), %r13	 #, tmp514
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	leaq	.LC10(%rip), %r14	 #, tmp511
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmulss	%xmm9, %xmm0, %xmm0	 # tmp519, tmp402, tmp404
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:226: 	    return _ToDur(static_cast<__to_rep>(
	vcvttss2siq	%xmm0, %rsi	 # tmp404, _120
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r12, %rcx	 # tmp510,
	movq	%rsi, %rdx	 # _120,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp524, _126
 # poly_eval_test.cpp:154:     std::cout << "default time = " << d.count() << "ms\n";
	movq	%r13, %rdx	 # tmp514,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$6, %r8d	 #,
	movq	%r14, %rdx	 # tmp511,
	movq	%r12, %rcx	 # tmp510,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_eval_test.cpp:155: 	std::cout << "sum = " << sum/Ntest <<std::endl;
	vmovss	.LC11(%rip), %xmm8	 #, tmp513
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp510,
 # poly_eval_test.cpp:155: 	std::cout << "sum = " << sum/Ntest <<std::endl;
	vdivss	%xmm8, %xmm12, %xmm1	 # tmp513, sum, tmp412
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # tmp412,
	call	_ZNSo9_M_insertIdEERSoT_	 #
	movq	%rax, %rcx	 # tmp525, _128
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0	 #
 # poly_eval_test.cpp:163: 	t0 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
	vmovss	_ZL9rand_init(%rip), %xmm0	 # rand_init, rand_init_lsm.463
	movq	%rax, %rdi	 # tmp526, t0
 # poly_eval_test.cpp:162: 	sum = 0;
	vmovaps	%xmm11, %xmm12	 # sum, sum
 # poly_eval_test.cpp:163: 	t0 = Time::now();
	movl	$100000000, %eax	 #, ivtmp_318
	.p2align 4,,10
	.p2align 3
.L18:
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm0, %xmm0, %xmm0	 # rand_init_lsm.463, tmp416
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm0, %xmm0	 # tmp516, tmp416, tmp417
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm0, %xmm0, %xmm0	 # tmp417, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm0, %xmm0, %xmm1	 #, y, _133
	vsubss	%xmm1, %xmm0, %xmm0	 # _133, y, rand_init_lsm.463
 # poly_eval_test.cpp:170: 		sum += rnd();
	vaddss	%xmm0, %xmm12, %xmm12	 # rand_init_lsm.463, sum, sum
 # poly_eval_test.cpp:164: 	for(size_t i=0;i<Ntest;++i){
	decq	%rax	 # ivtmp_318
	jne	.L18	 #,
	vmovss	%xmm0, _ZL9rand_init(%rip)	 # rand_init_lsm.463, rand_init
 # poly_eval_test.cpp:172: 	t1 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:653: 	return __cd(__cd(__lhs).count() - __cd(__rhs).count());
	subq	%rdi, %rax	 # t0, tmp419
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vcvtsi2ssq	%rax, %xmm7, %xmm0	 # tmp419, tmp548, tmp551
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$17, %r8d	 #,
	leaq	.LC12(%rip), %rdx	 #, tmp425
	movq	%r12, %rcx	 # tmp510,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vdivss	%xmm10, %xmm0, %xmm0	 # tmp518, tmp420, tmp421
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmulss	%xmm9, %xmm0, %xmm0	 # tmp519, tmp421, tmp423
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:226: 	    return _ToDur(static_cast<__to_rep>(
	vcvttss2siq	%xmm0, %rdi	 # tmp423, _136
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r12, %rcx	 # tmp510,
	movq	%rdi, %rdx	 # _136,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp528, _140
 # poly_eval_test.cpp:175: 	std::cout << "eval zero time = " << d.count() << "ms\n";
	movq	%r13, %rdx	 # tmp514,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$6, %r8d	 #,
	movq	%r14, %rdx	 # tmp511,
	movq	%r12, %rcx	 # tmp510,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_eval_test.cpp:176: 	std::cout << "sum = " << sum/Ntest <<std::endl;
	vdivss	%xmm8, %xmm12, %xmm1	 # tmp513, sum, tmp431
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp510,
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # tmp431,
	call	_ZNSo9_M_insertIdEERSoT_	 #
	movq	%rax, %rcx	 # tmp529, _142
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0	 #
 # poly_eval_test.cpp:181: 	t0 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # poly_eval_test.cpp:115: 	rand_init = next_f(rand_init);
	vmovss	_ZL9rand_init(%rip), %xmm1	 # rand_init, _148
 # poly_eval_test.cpp:181: 	t0 = Time::now();
	movq	%rax, %r15	 # tmp530, t0
 # poly_eval_test.cpp:180: 	sum = 0;
	vmovaps	%xmm11, %xmm12	 # sum, sum
 # poly_eval_test.cpp:115: 	rand_init = next_f(rand_init);
	movl	$100000000, %eax	 #, ivtmp_320
	.p2align 4,,10
	.p2align 3
.L19:
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # _148, tmp435
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm1, %xmm1	 # tmp516, tmp435, tmp436
 # poly_eval_test.cpp:188: 		sum += eval<8>(cf.data(),rnd());
	movq	%rbx, %rcx	 # tmp287,
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm1, %xmm1, %xmm1	 # tmp436, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm1, %xmm1, %xmm0	 #, y, _147
	vsubss	%xmm0, %xmm1, %xmm1	 # _147, y, _148
 # poly_eval_test.cpp:115: 	rand_init = next_f(rand_init);
	vmovss	%xmm1, _ZL9rand_init(%rip)	 # _148, rand_init
 # poly_eval_test.cpp:188: 		sum += eval<8>(cf.data(),rnd());
	vmovss	%xmm1, 32(%rsp)	 # _148, %sfp
	call	_Z4evalILy8EEfPff	 #
 # poly_eval_test.cpp:182: 	for(size_t i=0;i<Ntest;++i){
	decq	%rax	 # ivtmp_320
 # poly_eval_test.cpp:188: 		sum += eval<8>(cf.data(),rnd());
	vaddss	%xmm0, %xmm12, %xmm12	 # tmp531, sum, sum
 # poly_eval_test.cpp:182: 	for(size_t i=0;i<Ntest;++i){
	vmovss	32(%rsp), %xmm1	 # %sfp, _148
	jne	.L19	 #,
 # poly_eval_test.cpp:190: 	t1 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:653: 	return __cd(__cd(__lhs).count() - __cd(__rhs).count());
	subq	%r15, %rax	 # t0, tmp438
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vcvtsi2ssq	%rax, %xmm7, %xmm0	 # tmp438, tmp548, tmp552
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$21, %r8d	 #,
	leaq	.LC13(%rip), %rdx	 #, tmp444
	movq	%r12, %rcx	 # tmp510,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vdivss	%xmm10, %xmm0, %xmm0	 # tmp518, tmp439, tmp440
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	leaq	.LC14(%rip), %r15	 #, tmp512
	addq	$256, %rbx	 #, tmp515
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmulss	%xmm9, %xmm0, %xmm0	 # tmp519, tmp440, tmp442
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:226: 	    return _ToDur(static_cast<__to_rep>(
	vcvttss2siq	%xmm0, %rbp	 # tmp442, _151
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r12, %rcx	 # tmp510,
	movq	%rbp, %rdx	 # _151,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp533, _155
 # poly_eval_test.cpp:193: 	std::cout << "eval parallel time = " << d.count() << "ms\n";
	movq	%r13, %rdx	 # tmp514,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$13, %r8d	 #,
	movq	%r15, %rdx	 # tmp512,
	movq	%r12, %rcx	 # tmp510,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_eval_test.cpp:194: 	std::cout << "delta time = " << d.count() - def_time1  << "ms\n";
	movq	%rbp, %rdx	 # _151, _151
	subq	%rdi, %rdx	 # _136, _151
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r12, %rcx	 # tmp510,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp534, _156
 # poly_eval_test.cpp:194: 	std::cout << "delta time = " << d.count() - def_time1  << "ms\n";
	movq	%r13, %rdx	 # tmp514,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$6, %r8d	 #,
	movq	%r14, %rdx	 # tmp511,
	movq	%r12, %rcx	 # tmp510,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_eval_test.cpp:195: 	std::cout << "sum = " << sum/Ntest <<std::endl;
	vdivss	%xmm8, %xmm12, %xmm1	 # tmp513, sum, tmp455
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp510,
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # tmp455,
	call	_ZNSo9_M_insertIdEERSoT_	 #
	movq	%rax, %rcx	 # tmp535, _158
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0	 #
 # poly_eval_test.cpp:199: 	t0 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # poly_eval_test.cpp:115: 	rand_init = next_f(rand_init);
	vmovss	_ZL9rand_init(%rip), %xmm1	 # rand_init, _164
 # poly_eval_test.cpp:199: 	t0 = Time::now();
	movq	%rax, %rdi	 # tmp536, t0
 # poly_eval_test.cpp:198: 	sum = 0;
	vmovaps	%xmm11, %xmm12	 # sum, sum
 # poly_eval_test.cpp:115: 	rand_init = next_f(rand_init);
	movl	$100000000, %eax	 #, ivtmp_322
	.p2align 4,,10
	.p2align 3
.L20:
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # _164, tmp459
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm1, %xmm1	 # tmp516, tmp459, tmp460
 # poly_eval_test.cpp:206: 		sum += eval1<8>(coeffs,rnd());;
	movq	%rbx, %rcx	 # tmp515,
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm1, %xmm1, %xmm1	 # tmp460, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm1, %xmm1, %xmm0	 #, y, _163
	vsubss	%xmm0, %xmm1, %xmm1	 # _163, y, _164
 # poly_eval_test.cpp:115: 	rand_init = next_f(rand_init);
	vmovss	%xmm1, _ZL9rand_init(%rip)	 # _164, rand_init
 # poly_eval_test.cpp:206: 		sum += eval1<8>(coeffs,rnd());;
	call	_Z5eval1ILy8EEfPff	 #
 # poly_eval_test.cpp:206: 		sum += eval1<8>(coeffs,rnd());;
	vaddss	%xmm0, %xmm12, %xmm12	 # tmp537, sum, sum
 # poly_eval_test.cpp:200: 	for(size_t i=0;i<Ntest;++i){
	decq	%rax	 # ivtmp_322
	jne	.L20	 #,
 # poly_eval_test.cpp:208: 	t1 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:653: 	return __cd(__cd(__lhs).count() - __cd(__rhs).count());
	subq	%rdi, %rax	 # t0, tmp463
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vcvtsi2ssq	%rax, %xmm7, %xmm0	 # tmp463, tmp548, tmp553
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$25, %r8d	 #,
	leaq	.LC15(%rip), %rdx	 #, tmp469
	movq	%r12, %rcx	 # tmp510,
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vdivss	%xmm10, %xmm0, %xmm0	 # tmp518, tmp464, tmp465
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmulss	%xmm9, %xmm0, %xmm0	 # tmp519, tmp465, tmp467
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:226: 	    return _ToDur(static_cast<__to_rep>(
	vcvttss2siq	%xmm0, %rdi	 # tmp467, _166
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r12, %rcx	 # tmp510,
	movq	%rdi, %rdx	 # _166,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp539, _170
 # poly_eval_test.cpp:211: 	std::cout << "eval not parallel time = " << d.count() << "ms\n";
	movq	%r13, %rdx	 # tmp514,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0	 #
 # poly_eval_test.cpp:212: 	std::cout << "delta time = " << d.count() - def_time  << "ms\n";
	subq	%rsi, %rdi	 # _120, _166
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$13, %r8d	 #,
	movq	%r15, %rdx	 # tmp512,
	movq	%r12, %rcx	 # tmp510,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_eval_test.cpp:212: 	std::cout << "delta time = " << d.count() - def_time  << "ms\n";
	movq	%rdi, %rdx	 # _166, tmp475
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r12, %rcx	 # tmp510,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp540, _171
 # poly_eval_test.cpp:212: 	std::cout << "delta time = " << d.count() - def_time  << "ms\n";
	movq	%r13, %rdx	 # tmp514,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$6, %r8d	 #,
	movq	%r14, %rdx	 # tmp511,
	movq	%r12, %rcx	 # tmp510,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_eval_test.cpp:213: 	std::cout << "sum = " << sum/Ntest << std::endl;
	vdivss	%xmm8, %xmm12, %xmm1	 # tmp513, sum, tmp480
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp510,
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # tmp480,
	call	_ZNSo9_M_insertIdEERSoT_	 #
	movq	%rax, %rcx	 # tmp541, _173
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0	 #
 # poly_eval_test.cpp:216: 	t0 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # poly_eval_test.cpp:115: 	rand_init = next_f(rand_init);
	vmovss	_ZL9rand_init(%rip), %xmm1	 # rand_init, _179
 # poly_eval_test.cpp:216: 	t0 = Time::now();
	movq	%rax, %rdi	 # tmp542, t0
 # poly_eval_test.cpp:115: 	rand_init = next_f(rand_init);
	movl	$100000000, %eax	 #, ivtmp_324
	.p2align 4,,10
	.p2align 3
.L21:
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # _179, tmp484
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vaddsd	%xmm6, %xmm1, %xmm1	 # tmp516, tmp484, tmp485
 # poly_eval_test.cpp:223: 		sum += eval2<8>(coeffs,rnd());;
	movq	%rbx, %rcx	 # tmp515,
 # poly_eval_test.cpp:109: 	float y = x + 3.14;
	vcvtsd2ss	%xmm1, %xmm1, %xmm1	 # tmp485, y
 # poly_eval_test.cpp:110: 	return y - floor(y);
	vroundss	$9, %xmm1, %xmm1, %xmm0	 #, y, _178
	vsubss	%xmm0, %xmm1, %xmm1	 # _178, y, _179
 # poly_eval_test.cpp:115: 	rand_init = next_f(rand_init);
	vmovss	%xmm1, _ZL9rand_init(%rip)	 # _179, rand_init
 # poly_eval_test.cpp:223: 		sum += eval2<8>(coeffs,rnd());;
	call	_Z5eval2ILy8EEfPff	 #
 # poly_eval_test.cpp:223: 		sum += eval2<8>(coeffs,rnd());;
	vaddss	%xmm0, %xmm11, %xmm11	 # tmp543, sum, sum
 # poly_eval_test.cpp:217: 	for(size_t i=0;i<Ntest;++i){
	decq	%rax	 # ivtmp_324
	jne	.L21	 #,
 # poly_eval_test.cpp:225: 	t1 = Time::now();
	call	_ZNSt6chrono3_V212system_clock3nowEv	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:653: 	return __cd(__cd(__lhs).count() - __cd(__rhs).count());
	subq	%rdi, %rax	 # t0, tmp488
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vcvtsi2ssq	%rax, %xmm7, %xmm7	 # tmp488, tmp548, tmp554
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$25, %r8d	 #,
	movq	%r12, %rcx	 # tmp510,
	leaq	.LC16(%rip), %rdx	 #, tmp494
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:214: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	vdivss	%xmm10, %xmm7, %xmm0	 # tmp518, tmp489, tmp490
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:227: 	      static_cast<_CR>(__d.count()) * static_cast<_CR>(_CF::num)));
	vmulss	%xmm9, %xmm0, %xmm0	 # tmp519, tmp490, tmp492
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/chrono:226: 	    return _ToDur(static_cast<__to_rep>(
	vcvttss2siq	%xmm0, %rbx	 # tmp492, _181
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r12, %rcx	 # tmp510,
	movq	%rbx, %rdx	 # _181,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp545, _185
 # poly_eval_test.cpp:228: 	std::cout << "not simd parallel time = " << d.count() << "ms\n";
	movq	%r13, %rdx	 # tmp514,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$13, %r8d	 #,
	movq	%r15, %rdx	 # tmp512,
	movq	%r12, %rcx	 # tmp510,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_eval_test.cpp:229: 	std::cout << "delta time = " << d.count() - def_time  << "ms\n";
	movq	%rbx, %rdx	 # _181, _181
	subq	%rsi, %rdx	 # _120, _181
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:202:       { return _M_insert(__n); }
	movq	%r12, %rcx	 # tmp510,
	call	_ZNSo9_M_insertIxEERSoT_	 #
	movq	%rax, %rcx	 # tmp546, _186
 # poly_eval_test.cpp:229: 	std::cout << "delta time = " << d.count() - def_time  << "ms\n";
	movq	%r13, %rdx	 # tmp514,
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0	 #
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:611: 	__ostream_insert(__out, __s,
	movl	$6, %r8d	 #,
	movq	%r14, %rdx	 # tmp511,
	movq	%r12, %rcx	 # tmp510,
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x	 #
 # poly_eval_test.cpp:230: 	std::cout << "sum = " << sum/Ntest << std::endl;
	vdivss	%xmm8, %xmm11, %xmm1	 # tmp513, sum, tmp505
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:228: 	return _M_insert(static_cast<double>(__f));
	movq	%r12, %rcx	 # tmp510,
	vcvtss2sd	%xmm1, %xmm1, %xmm1	 # tmp505,
	call	_ZNSo9_M_insertIdEERSoT_	 #
	movq	%rax, %rcx	 # tmp547, _188
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/ostream:113: 	return __pf(*this);
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0	 #
	nop	
 # poly_eval_test.cpp:232: }
	vmovaps	592(%rsp), %xmm6	 #,
	vmovaps	608(%rsp), %xmm7	 #,
	vmovaps	624(%rsp), %xmm8	 #,
	vmovaps	640(%rsp), %xmm9	 #,
	vmovaps	656(%rsp), %xmm10	 #,
	vmovaps	672(%rsp), %xmm11	 #,
	vmovaps	688(%rsp), %xmm12	 #,
	vmovaps	704(%rsp), %xmm13	 #,
	vmovaps	720(%rsp), %xmm14	 #,
	vmovaps	736(%rsp), %xmm15	 #,
	xorl	%eax, %eax	 #
	addq	$760, %rsp	 #,
	popq	%rbx	 #
	popq	%rsi	 #
	popq	%rdi	 #
	popq	%rbp	 #
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
.LFB14746:
	subq	$40, %rsp	 #,
	.seh_stackalloc	40
	.seh_endprologue
 # D:/Soft/Qt/Tools/mingw1120_64/lib/gcc/x86_64-w64-mingw32/11.2.0/include/c++/iostream:74:   static ios_base::Init __ioinit;
	leaq	_ZStL8__ioinit(%rip), %rcx	 #, tmp82
	call	_ZNSt8ios_base4InitC1Ev	 #
	leaq	__tcf_1(%rip), %rcx	 #, tmp83
 # poly_eval_test.cpp:232: }
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
	.set	.LC0,.LC5
	.section .rdata,"dr"
	.align 32
.LC2:
	.long	1065353216
	.long	1073741824
	.long	1077936128
	.long	1082130432
	.long	1084227584
	.long	1086324736
	.long	1088421888
	.long	1090519040
	.align 32
.LC3:
	.long	1065353216
	.long	1084227584
	.long	1077936128
	.long	1088421888
	.long	1073741824
	.long	1086324736
	.long	1082130432
	.long	1090519040
	.align 4
.LC4:
	.long	1191181824
	.align 32
.LC5:
	.long	1374389535
	.long	1074339512
	.long	1374389535
	.long	1074339512
	.long	1374389535
	.long	1074339512
	.long	1374389535
	.long	1074339512
	.align 4
.LC6:
	.long	1315859240
	.align 4
.LC7:
	.long	1148846080
	.align 4
.LC11:
	.long	1287568416
	.ident	"GCC: (x86_64-posix-seh-rev3, Built by MinGW-W64 project) 11.2.0"
	.def	_ZNSt8ios_base4InitD1Ev;	.scl	2;	.type	32;	.endef
	.def	_ZNSo3putEc;	.scl	2;	.type	32;	.endef
	.def	_ZNSo5flushEv;	.scl	2;	.type	32;	.endef
	.def	_ZNKSt5ctypeIcE13_M_widen_initEv;	.scl	2;	.type	32;	.endef
	.def	_ZSt16__throw_bad_castv;	.scl	2;	.type	32;	.endef
	.def	strlen;	.scl	2;	.type	32;	.endef
	.def	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_x;	.scl	2;	.type	32;	.endef
	.def	_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate;	.scl	2;	.type	32;	.endef
	.def	srand;	.scl	2;	.type	32;	.endef
	.def	rand;	.scl	2;	.type	32;	.endef
	.def	_ZNSt6chrono3_V212system_clock3nowEv;	.scl	2;	.type	32;	.endef
	.def	_ZNSo9_M_insertIxEERSoT_;	.scl	2;	.type	32;	.endef
	.def	_ZNSo9_M_insertIdEERSoT_;	.scl	2;	.type	32;	.endef
	.def	_ZNSt8ios_base4InitC1Ev;	.scl	2;	.type	32;	.endef
	.def	atexit;	.scl	2;	.type	32;	.endef
	.section	.rdata$.refptr._ZSt4cout, "dr"
	.globl	.refptr._ZSt4cout
	.linkonce	discard
.refptr._ZSt4cout:
	.quad	_ZSt4cout
