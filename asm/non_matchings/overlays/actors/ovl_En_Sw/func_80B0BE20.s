glabel func_80B0BE20
/* 00070 80B0BE20 27BDFFB0 */  addiu   $sp, $sp, 0xFFB0           ## $sp = FFFFFFB0
/* 00074 80B0BE24 3C0180B1 */  lui     $at, %hi(D_80B0F200)       ## $at = 80B10000
/* 00078 80B0BE28 C420F200 */  lwc1    $f0, %lo(D_80B0F200)($at)  
/* 0007C 80B0BE2C AFBF001C */  sw      $ra, 0x001C($sp)           
/* 00080 80B0BE30 AFB00018 */  sw      $s0, 0x0018($sp)           
/* 00084 80B0BE34 AC850078 */  sw      $a1, 0x0078($a0)           ## 00000078
/* 00088 80B0BE38 84AE0008 */  lh      $t6, 0x0008($a1)           ## 00000008
/* 0008C 80B0BE3C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00090 80B0BE40 448E2000 */  mtc1    $t6, $f4                   ## $f4 = 0.00
/* 00094 80B0BE44 00000000 */  nop
/* 00098 80B0BE48 468021A0 */  cvt.s.w $f6, $f4                   
/* 0009C 80B0BE4C 46003202 */  mul.s   $f8, $f6, $f0              
/* 000A0 80B0BE50 E7A80044 */  swc1    $f8, 0x0044($sp)           
/* 000A4 80B0BE54 84AF000A */  lh      $t7, 0x000A($a1)           ## 0000000A
/* 000A8 80B0BE58 448F5000 */  mtc1    $t7, $f10                  ## $f10 = 0.00
/* 000AC 80B0BE5C 00000000 */  nop
/* 000B0 80B0BE60 46805420 */  cvt.s.w $f16, $f10                 
/* 000B4 80B0BE64 46008482 */  mul.s   $f18, $f16, $f0            
/* 000B8 80B0BE68 E7B20048 */  swc1    $f18, 0x0048($sp)          
/* 000BC 80B0BE6C 84B8000C */  lh      $t8, 0x000C($a1)           ## 0000000C
/* 000C0 80B0BE70 44982000 */  mtc1    $t8, $f4                   ## $f4 = 0.00
/* 000C4 80B0BE74 00000000 */  nop
/* 000C8 80B0BE78 468021A0 */  cvt.s.w $f6, $f4                   
/* 000CC 80B0BE7C 46003282 */  mul.s   $f10, $f6, $f0             
/* 000D0 80B0BE80 E7AA004C */  swc1    $f10, 0x004C($sp)          
/* 000D4 80B0BE84 C4900364 */  lwc1    $f16, 0x0364($a0)          ## 00000364
/* 000D8 80B0BE88 C4860368 */  lwc1    $f6, 0x0368($a0)           ## 00000368
/* 000DC 80B0BE8C 46088102 */  mul.s   $f4, $f16, $f8             
/* 000E0 80B0BE90 00000000 */  nop
/* 000E4 80B0BE94 46069402 */  mul.s   $f16, $f18, $f6            
/* 000E8 80B0BE98 C492036C */  lwc1    $f18, 0x036C($a0)          ## 0000036C
/* 000EC 80B0BE9C 46125182 */  mul.s   $f6, $f10, $f18            
/* 000F0 80B0BEA0 46102200 */  add.s   $f8, $f4, $f16             
/* 000F4 80B0BEA4 0C03F4DA */  jal     func_800FD368              
/* 000F8 80B0BEA8 46064300 */  add.s   $f12, $f8, $f6             
/* 000FC 80B0BEAC 26040364 */  addiu   $a0, $s0, 0x0364           ## $a0 = 00000364
/* 00100 80B0BEB0 E7A00034 */  swc1    $f0, 0x0034($sp)           
/* 00104 80B0BEB4 AFA40028 */  sw      $a0, 0x0028($sp)           
/* 00108 80B0BEB8 27A50044 */  addiu   $a1, $sp, 0x0044           ## $a1 = FFFFFFF4
/* 0010C 80B0BEBC 0C2C2F6C */  jal     func_80B0BDB0              
/* 00110 80B0BEC0 27A60038 */  addiu   $a2, $sp, 0x0038           ## $a2 = FFFFFFE8
/* 00114 80B0BEC4 C7AC0034 */  lwc1    $f12, 0x0034($sp)          
/* 00118 80B0BEC8 27A50038 */  addiu   $a1, $sp, 0x0038           ## $a1 = FFFFFFE8
/* 0011C 80B0BECC 0C0348FF */  jal     func_800D23FC              
/* 00120 80B0BED0 00003025 */  or      $a2, $zero, $zero          ## $a2 = 00000000
/* 00124 80B0BED4 26040370 */  addiu   $a0, $s0, 0x0370           ## $a0 = 00000370
/* 00128 80B0BED8 AFA40024 */  sw      $a0, 0x0024($sp)           
/* 0012C 80B0BEDC 0C0346BD */  jal     Matrix_MultVec3f              
/* 00130 80B0BEE0 27A50038 */  addiu   $a1, $sp, 0x0038           ## $a1 = FFFFFFE8
/* 00134 80B0BEE4 27B90038 */  addiu   $t9, $sp, 0x0038           ## $t9 = FFFFFFE8
/* 00138 80B0BEE8 8FA40024 */  lw      $a0, 0x0024($sp)           
/* 0013C 80B0BEEC 8F290000 */  lw      $t1, 0x0000($t9)           ## FFFFFFE8
/* 00140 80B0BEF0 2606037C */  addiu   $a2, $s0, 0x037C           ## $a2 = 0000037C
/* 00144 80B0BEF4 27A50044 */  addiu   $a1, $sp, 0x0044           ## $a1 = FFFFFFF4
/* 00148 80B0BEF8 AC890000 */  sw      $t1, 0x0000($a0)           ## 00000000
/* 0014C 80B0BEFC 8F280004 */  lw      $t0, 0x0004($t9)           ## FFFFFFEC
/* 00150 80B0BF00 AC880004 */  sw      $t0, 0x0004($a0)           ## 00000004
/* 00154 80B0BF04 8F290008 */  lw      $t1, 0x0008($t9)           ## FFFFFFF0
/* 00158 80B0BF08 AC890008 */  sw      $t1, 0x0008($a0)           ## 00000008
/* 0015C 80B0BF0C 0C2C2F6C */  jal     func_80B0BDB0              
/* 00160 80B0BF10 AFA60020 */  sw      $a2, 0x0020($sp)           
/* 00164 80B0BF14 0C032D8A */  jal     func_800CB628              
/* 00168 80B0BF18 8FA40020 */  lw      $a0, 0x0020($sp)           
/* 0016C 80B0BF1C 3C0180B1 */  lui     $at, %hi(D_80B0F204)       ## $at = 80B10000
/* 00170 80B0BF20 C424F204 */  lwc1    $f4, %lo(D_80B0F204)($at)  
/* 00174 80B0BF24 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 00178 80B0BF28 4604003C */  c.lt.s  $f0, $f4                   
/* 0017C 80B0BF2C 00000000 */  nop
/* 00180 80B0BF30 45020004 */  bc1fl   .L80B0BF44                 
/* 00184 80B0BF34 44817000 */  mtc1    $at, $f14                  ## $f14 = 1.00
/* 00188 80B0BF38 10000034 */  beq     $zero, $zero, .L80B0C00C   
/* 0018C 80B0BF3C 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 00190 80B0BF40 44817000 */  mtc1    $at, $f14                  ## $f14 = 1.00
.L80B0BF44:
/* 00194 80B0BF44 C610037C */  lwc1    $f16, 0x037C($s0)          ## 0000037C
/* 00198 80B0BF48 C6120380 */  lwc1    $f18, 0x0380($s0)          ## 00000380
/* 0019C 80B0BF4C 46007083 */  div.s   $f2, $f14, $f0             
/* 001A0 80B0BF50 C6060384 */  lwc1    $f6, 0x0384($s0)           ## 00000384
/* 001A4 80B0BF54 27AB0044 */  addiu   $t3, $sp, 0x0044           ## $t3 = FFFFFFF4
/* 001A8 80B0BF58 44806000 */  mtc1    $zero, $f12                ## $f12 = 0.00
/* 001AC 80B0BF5C 260403D8 */  addiu   $a0, $s0, 0x03D8           ## $a0 = 000003D8
/* 001B0 80B0BF60 26050030 */  addiu   $a1, $s0, 0x0030           ## $a1 = 00000030
/* 001B4 80B0BF64 00003025 */  or      $a2, $zero, $zero          ## $a2 = 00000000
/* 001B8 80B0BF68 46028282 */  mul.s   $f10, $f16, $f2            
/* 001BC 80B0BF6C 00000000 */  nop
/* 001C0 80B0BF70 46029202 */  mul.s   $f8, $f18, $f2             
/* 001C4 80B0BF74 00000000 */  nop
/* 001C8 80B0BF78 46023102 */  mul.s   $f4, $f6, $f2              
/* 001CC 80B0BF7C E60A037C */  swc1    $f10, 0x037C($s0)          ## 0000037C
/* 001D0 80B0BF80 E6080380 */  swc1    $f8, 0x0380($s0)           ## 00000380
/* 001D4 80B0BF84 E6040384 */  swc1    $f4, 0x0384($s0)           ## 00000384
/* 001D8 80B0BF88 8D6D0000 */  lw      $t5, 0x0000($t3)           ## FFFFFFF4
/* 001DC 80B0BF8C 8FAA0028 */  lw      $t2, 0x0028($sp)           
/* 001E0 80B0BF90 AD4D0000 */  sw      $t5, 0x0000($t2)           ## 00000000
/* 001E4 80B0BF94 8D6C0004 */  lw      $t4, 0x0004($t3)           ## FFFFFFF8
/* 001E8 80B0BF98 AD4C0004 */  sw      $t4, 0x0004($t2)           ## 00000004
/* 001EC 80B0BF9C 8D6D0008 */  lw      $t5, 0x0008($t3)           ## FFFFFFFC
/* 001F0 80B0BFA0 AD4D0008 */  sw      $t5, 0x0008($t2)           ## 00000008
/* 001F4 80B0BFA4 C6100370 */  lwc1    $f16, 0x0370($s0)          ## 00000370
/* 001F8 80B0BFA8 C60A0374 */  lwc1    $f10, 0x0374($s0)          ## 00000374
/* 001FC 80B0BFAC C6120378 */  lwc1    $f18, 0x0378($s0)          ## 00000378
/* 00200 80B0BFB0 E61003D8 */  swc1    $f16, 0x03D8($s0)          ## 000003D8
/* 00204 80B0BFB4 E60A03DC */  swc1    $f10, 0x03DC($s0)          ## 000003DC
/* 00208 80B0BFB8 E61203E0 */  swc1    $f18, 0x03E0($s0)          ## 000003E0
/* 0020C 80B0BFBC C6080364 */  lwc1    $f8, 0x0364($s0)           ## 00000364
/* 00210 80B0BFC0 C6060368 */  lwc1    $f6, 0x0368($s0)           ## 00000368
/* 00214 80B0BFC4 C604036C */  lwc1    $f4, 0x036C($s0)           ## 0000036C
/* 00218 80B0BFC8 C610037C */  lwc1    $f16, 0x037C($s0)          ## 0000037C
/* 0021C 80B0BFCC C60A0380 */  lwc1    $f10, 0x0380($s0)          ## 00000380
/* 00220 80B0BFD0 C6120384 */  lwc1    $f18, 0x0384($s0)          ## 00000384
/* 00224 80B0BFD4 E60E0414 */  swc1    $f14, 0x0414($s0)          ## 00000414
/* 00228 80B0BFD8 E60C03E4 */  swc1    $f12, 0x03E4($s0)          ## 000003E4
/* 0022C 80B0BFDC E60C03F4 */  swc1    $f12, 0x03F4($s0)          ## 000003F4
/* 00230 80B0BFE0 E60C0404 */  swc1    $f12, 0x0404($s0)          ## 00000404
/* 00234 80B0BFE4 E60C0408 */  swc1    $f12, 0x0408($s0)          ## 00000408
/* 00238 80B0BFE8 E60C040C */  swc1    $f12, 0x040C($s0)          ## 0000040C
/* 0023C 80B0BFEC E60C0410 */  swc1    $f12, 0x0410($s0)          ## 00000410
/* 00240 80B0BFF0 E60803E8 */  swc1    $f8, 0x03E8($s0)           ## 000003E8
/* 00244 80B0BFF4 E60603EC */  swc1    $f6, 0x03EC($s0)           ## 000003EC
/* 00248 80B0BFF8 E60403F0 */  swc1    $f4, 0x03F0($s0)           ## 000003F0
/* 0024C 80B0BFFC E61003F8 */  swc1    $f16, 0x03F8($s0)          ## 000003F8
/* 00250 80B0C000 E60A03FC */  swc1    $f10, 0x03FC($s0)          ## 000003FC
/* 00254 80B0C004 0C034833 */  jal     func_800D20CC              
/* 00258 80B0C008 E6120400 */  swc1    $f18, 0x0400($s0)          ## 00000400
.L80B0C00C:
/* 0025C 80B0C00C 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 00260 80B0C010 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 00264 80B0C014 27BD0050 */  addiu   $sp, $sp, 0x0050           ## $sp = 00000000
/* 00268 80B0C018 03E00008 */  jr      $ra                        
/* 0026C 80B0C01C 00000000 */  nop


