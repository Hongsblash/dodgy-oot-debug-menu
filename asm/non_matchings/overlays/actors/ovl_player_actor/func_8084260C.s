glabel func_8084260C
/* 103FC 8084260C 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 10400 80842610 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 10404 80842614 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 10408 80842618 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 1040C 8084261C AFA60020 */  sw      $a2, 0x0020($sp)           
/* 10410 80842620 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 10414 80842624 AFA70024 */  sw      $a3, 0x0024($sp)           
/* 10418 80842628 C7A40024 */  lwc1    $f4, 0x0024($sp)           
/* 1041C 8084262C 8FAE0018 */  lw      $t6, 0x0018($sp)           
/* 10420 80842630 8FAF001C */  lw      $t7, 0x001C($sp)           
/* 10424 80842634 46040182 */  mul.s   $f6, $f0, $f4              
/* 10428 80842638 C5C80000 */  lwc1    $f8, 0x0000($t6)           ## 00000000
/* 1042C 8084263C 46083280 */  add.s   $f10, $f6, $f8             
/* 10430 80842640 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 10434 80842644 E5EA0000 */  swc1    $f10, 0x0000($t7)          ## 00000000
/* 10438 80842648 8FB80018 */  lw      $t8, 0x0018($sp)           
/* 1043C 8084264C C7A60028 */  lwc1    $f6, 0x0028($sp)           
/* 10440 80842650 C7B20020 */  lwc1    $f18, 0x0020($sp)          
/* 10444 80842654 C7100004 */  lwc1    $f16, 0x0004($t8)          ## 00000004
/* 10448 80842658 46060202 */  mul.s   $f8, $f0, $f6              
/* 1044C 8084265C 8FB9001C */  lw      $t9, 0x001C($sp)           
/* 10450 80842660 46128100 */  add.s   $f4, $f16, $f18            
/* 10454 80842664 46044280 */  add.s   $f10, $f8, $f4             
/* 10458 80842668 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 1045C 8084266C E72A0004 */  swc1    $f10, 0x0004($t9)          ## 00000004
/* 10460 80842670 C7B00024 */  lwc1    $f16, 0x0024($sp)          
/* 10464 80842674 8FA80018 */  lw      $t0, 0x0018($sp)           
/* 10468 80842678 8FA9001C */  lw      $t1, 0x001C($sp)           
/* 1046C 8084267C 46100482 */  mul.s   $f18, $f0, $f16            
/* 10470 80842680 C5060008 */  lwc1    $f6, 0x0008($t0)           ## 00000008
/* 10474 80842684 46069200 */  add.s   $f8, $f18, $f6             
/* 10478 80842688 E5280008 */  swc1    $f8, 0x0008($t1)           ## 00000008
/* 1047C 8084268C 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 10480 80842690 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 10484 80842694 03E00008 */  jr      $ra                        
/* 10488 80842698 00000000 */  nop
