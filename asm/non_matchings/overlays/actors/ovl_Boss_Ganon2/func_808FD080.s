glabel func_808FD080
/* 00140 808FD080 C4C40000 */  lwc1    $f4, 0x0000($a2)           ## 00000000
/* 00144 808FD084 8CB8001C */  lw      $t8, 0x001C($a1)           ## 0000001C
/* 00148 808FD088 00041180 */  sll     $v0, $a0,  6               
/* 0014C 808FD08C 4600218D */  trunc.w.s $f6, $f4                   
/* 00150 808FD090 0302C821 */  addu    $t9, $t8, $v0              
/* 00154 808FD094 440F3000 */  mfc1    $t7, $f6                   
/* 00158 808FD098 00000000 */  nop
/* 0015C 808FD09C A72F0030 */  sh      $t7, 0x0030($t9)           ## 00000030
/* 00160 808FD0A0 C4C80004 */  lwc1    $f8, 0x0004($a2)           ## 00000004
/* 00164 808FD0A4 8CAA001C */  lw      $t2, 0x001C($a1)           ## 0000001C
/* 00168 808FD0A8 4600428D */  trunc.w.s $f10, $f8                  
/* 0016C 808FD0AC 01425821 */  addu    $t3, $t2, $v0              
/* 00170 808FD0B0 44095000 */  mfc1    $t1, $f10                  
/* 00174 808FD0B4 00000000 */  nop
/* 00178 808FD0B8 A5690032 */  sh      $t1, 0x0032($t3)           ## 00000032
/* 0017C 808FD0BC C4D00008 */  lwc1    $f16, 0x0008($a2)          ## 00000008
/* 00180 808FD0C0 8CAE001C */  lw      $t6, 0x001C($a1)           ## 0000001C
/* 00184 808FD0C4 4600848D */  trunc.w.s $f18, $f16                 
/* 00188 808FD0C8 01C2C021 */  addu    $t8, $t6, $v0              
/* 0018C 808FD0CC 440D9000 */  mfc1    $t5, $f18                  
/* 00190 808FD0D0 00000000 */  nop
/* 00194 808FD0D4 A70D0034 */  sh      $t5, 0x0034($t8)           ## 00000034
/* 00198 808FD0D8 8CAF001C */  lw      $t7, 0x001C($a1)           ## 0000001C
/* 0019C 808FD0DC 01E21821 */  addu    $v1, $t7, $v0              
/* 001A0 808FD0E0 8479002E */  lh      $t9, 0x002E($v1)           ## 0000002E
/* 001A4 808FD0E4 C4640038 */  lwc1    $f4, 0x0038($v1)           ## 00000038
/* 001A8 808FD0E8 44993000 */  mtc1    $t9, $f6                   ## $f6 = 0.00
/* 001AC 808FD0EC 00000000 */  nop
/* 001B0 808FD0F0 46803220 */  cvt.s.w $f8, $f6                   
/* 001B4 808FD0F4 46082282 */  mul.s   $f10, $f4, $f8             
/* 001B8 808FD0F8 4600540D */  trunc.w.s $f16, $f10                 
/* 001BC 808FD0FC 440A8000 */  mfc1    $t2, $f16                  
/* 001C0 808FD100 03E00008 */  jr      $ra                        
/* 001C4 808FD104 A46A0036 */  sh      $t2, 0x0036($v1)           ## 00000036


