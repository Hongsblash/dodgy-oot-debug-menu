glabel func_80AB08A4
/* 00484 80AB08A4 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 00488 80AB08A8 24010005 */  addiu   $at, $zero, 0x0005         ## $at = 00000005
/* 0048C 80AB08AC AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00490 80AB08B0 AFA40028 */  sw      $a0, 0x0028($sp)           
/* 00494 80AB08B4 AFA60030 */  sw      $a2, 0x0030($sp)           
/* 00498 80AB08B8 10A1000E */  beq     $a1, $at, .L80AB08F4       
/* 0049C 80AB08BC AFA70034 */  sw      $a3, 0x0034($sp)           
/* 004A0 80AB08C0 24010006 */  addiu   $at, $zero, 0x0006         ## $at = 00000006
/* 004A4 80AB08C4 10A1000B */  beq     $a1, $at, .L80AB08F4       
/* 004A8 80AB08C8 24010007 */  addiu   $at, $zero, 0x0007         ## $at = 00000007
/* 004AC 80AB08CC 10A10009 */  beq     $a1, $at, .L80AB08F4       
/* 004B0 80AB08D0 2401000B */  addiu   $at, $zero, 0x000B         ## $at = 0000000B
/* 004B4 80AB08D4 10A10007 */  beq     $a1, $at, .L80AB08F4       
/* 004B8 80AB08D8 2401000C */  addiu   $at, $zero, 0x000C         ## $at = 0000000C
/* 004BC 80AB08DC 10A10005 */  beq     $a1, $at, .L80AB08F4       
/* 004C0 80AB08E0 2401000D */  addiu   $at, $zero, 0x000D         ## $at = 0000000D
/* 004C4 80AB08E4 10A10003 */  beq     $a1, $at, .L80AB08F4       
/* 004C8 80AB08E8 2401000E */  addiu   $at, $zero, 0x000E         ## $at = 0000000E
/* 004CC 80AB08EC 54A10025 */  bnel    $a1, $at, .L80AB0984       
/* 004D0 80AB08F0 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80AB08F4:
/* 004D4 80AB08F4 8FAE003C */  lw      $t6, 0x003C($sp)           
/* 004D8 80AB08F8 00057840 */  sll     $t7, $a1,  1               
/* 004DC 80AB08FC 01CF1821 */  addu    $v1, $t6, $t7              
/* 004E0 80AB0900 8464020A */  lh      $a0, 0x020A($v1)           ## 0000020A
/* 004E4 80AB0904 0C01DE1C */  jal     Math_Sins
              ## sins?
/* 004E8 80AB0908 AFA3001C */  sw      $v1, 0x001C($sp)           
/* 004EC 80AB090C 8FA20038 */  lw      $v0, 0x0038($sp)           
/* 004F0 80AB0910 3C014348 */  lui     $at, 0x4348                ## $at = 43480000
/* 004F4 80AB0914 44814000 */  mtc1    $at, $f8                   ## $f8 = 200.00
/* 004F8 80AB0918 84580002 */  lh      $t8, 0x0002($v0)           ## 00000002
/* 004FC 80AB091C 8FA3001C */  lw      $v1, 0x001C($sp)           
/* 00500 80AB0920 46080282 */  mul.s   $f10, $f0, $f8             
/* 00504 80AB0924 44982000 */  mtc1    $t8, $f4                   ## $f4 = 0.00
/* 00508 80AB0928 00000000 */  nop
/* 0050C 80AB092C 468021A0 */  cvt.s.w $f6, $f4                   
/* 00510 80AB0930 460A3400 */  add.s   $f16, $f6, $f10            
/* 00514 80AB0934 4600848D */  trunc.w.s $f18, $f16                 
/* 00518 80AB0938 44089000 */  mfc1    $t0, $f18                  
/* 0051C 80AB093C 00000000 */  nop
/* 00520 80AB0940 A4480002 */  sh      $t0, 0x0002($v0)           ## 00000002
/* 00524 80AB0944 0C01DE0D */  jal     Math_Coss
              ## coss?
/* 00528 80AB0948 8464022A */  lh      $a0, 0x022A($v1)           ## 0000022A
/* 0052C 80AB094C 8FA20038 */  lw      $v0, 0x0038($sp)           
/* 00530 80AB0950 3C014348 */  lui     $at, 0x4348                ## $at = 43480000
/* 00534 80AB0954 44813000 */  mtc1    $at, $f6                   ## $f6 = 200.00
/* 00538 80AB0958 84490004 */  lh      $t1, 0x0004($v0)           ## 00000004
/* 0053C 80AB095C 46060282 */  mul.s   $f10, $f0, $f6             
/* 00540 80AB0960 44892000 */  mtc1    $t1, $f4                   ## $f4 = 0.00
/* 00544 80AB0964 00000000 */  nop
/* 00548 80AB0968 46802220 */  cvt.s.w $f8, $f4                   
/* 0054C 80AB096C 460A4400 */  add.s   $f16, $f8, $f10            
/* 00550 80AB0970 4600848D */  trunc.w.s $f18, $f16                 
/* 00554 80AB0974 440B9000 */  mfc1    $t3, $f18                  
/* 00558 80AB0978 00000000 */  nop
/* 0055C 80AB097C A44B0004 */  sh      $t3, 0x0004($v0)           ## 00000004
/* 00560 80AB0980 8FBF0014 */  lw      $ra, 0x0014($sp)           
.L80AB0984:
/* 00564 80AB0984 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 00568 80AB0988 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 0056C 80AB098C 03E00008 */  jr      $ra                        
/* 00570 80AB0990 00000000 */  nop


