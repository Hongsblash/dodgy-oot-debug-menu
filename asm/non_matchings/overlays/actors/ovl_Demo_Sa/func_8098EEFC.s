glabel func_8098EEFC
/* 00A9C 8098EEFC 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00AA0 8098EF00 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00AA4 8098EF04 AFA40018 */  sw      $a0, 0x0018($sp)           
/* 00AA8 8098EF08 24060004 */  addiu   $a2, $zero, 0x0004         ## $a2 = 00000004
/* 00AAC 8098EF0C 0C263995 */  jal     func_8098E654              
/* 00AB0 8098EF10 24070004 */  addiu   $a3, $zero, 0x0004         ## $a3 = 00000004
/* 00AB4 8098EF14 10400025 */  beq     $v0, $zero, .L8098EFAC     
/* 00AB8 8098EF18 8FA40018 */  lw      $a0, 0x0018($sp)           
/* 00ABC 8098EF1C 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 00AC0 8098EF20 44811000 */  mtc1    $at, $f2                   ## $f2 = 10.00
/* 00AC4 8098EF24 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 00AC8 8098EF28 248201A0 */  addiu   $v0, $a0, 0x01A0           ## $v0 = 000001A0
/* 00ACC 8098EF2C C4440000 */  lwc1    $f4, 0x0000($v0)           ## 000001A0
/* 00AD0 8098EF30 44813000 */  mtc1    $at, $f6                   ## $f6 = 1.00
/* 00AD4 8098EF34 3C058016 */  lui     $a1, 0x8016                ## $a1 = 80160000
/* 00AD8 8098EF38 24A5FA90 */  addiu   $a1, $a1, 0xFA90           ## $a1 = 8015FA90
/* 00ADC 8098EF3C 46062200 */  add.s   $f8, $f4, $f6              
/* 00AE0 8098EF40 24180009 */  addiu   $t8, $zero, 0x0009         ## $t8 = 00000009
/* 00AE4 8098EF44 24190001 */  addiu   $t9, $zero, 0x0001         ## $t9 = 00000001
/* 00AE8 8098EF48 E4480000 */  swc1    $f8, 0x0000($v0)           ## 000001A0
/* 00AEC 8098EF4C 8CAE0000 */  lw      $t6, 0x0000($a1)           ## 8015FA90
/* 00AF0 8098EF50 C4400000 */  lwc1    $f0, 0x0000($v0)           ## 000001A0
/* 00AF4 8098EF54 85CF145E */  lh      $t7, 0x145E($t6)           ## 0000145E
/* 00AF8 8098EF58 448F5000 */  mtc1    $t7, $f10                  ## $f10 = 0.00
/* 00AFC 8098EF5C 00000000 */  nop
/* 00B00 8098EF60 46805420 */  cvt.s.w $f16, $f10                 
/* 00B04 8098EF64 46028480 */  add.s   $f18, $f16, $f2            
/* 00B08 8098EF68 4600903E */  c.le.s  $f18, $f0                  
/* 00B0C 8098EF6C 00000000 */  nop
/* 00B10 8098EF70 45000020 */  bc1f    .L8098EFF4                 
/* 00B14 8098EF74 00000000 */  nop
/* 00B18 8098EF78 AC980198 */  sw      $t8, 0x0198($a0)           ## 00000198
/* 00B1C 8098EF7C AC99019C */  sw      $t9, 0x019C($a0)           ## 0000019C
/* 00B20 8098EF80 8CA80000 */  lw      $t0, 0x0000($a1)           ## 8015FA90
/* 00B24 8098EF84 240300FF */  addiu   $v1, $zero, 0x00FF         ## $v1 = 000000FF
/* 00B28 8098EF88 8509145E */  lh      $t1, 0x145E($t0)           ## 0000145E
/* 00B2C 8098EF8C 44892000 */  mtc1    $t1, $f4                   ## $f4 = 0.00
/* 00B30 8098EF90 00000000 */  nop
/* 00B34 8098EF94 468021A0 */  cvt.s.w $f6, $f4                   
/* 00B38 8098EF98 46023200 */  add.s   $f8, $f6, $f2              
/* 00B3C 8098EF9C E4480000 */  swc1    $f8, 0x0000($v0)           ## 000001A0
/* 00B40 8098EFA0 AC8301A4 */  sw      $v1, 0x01A4($a0)           ## 000001A4
/* 00B44 8098EFA4 10000026 */  beq     $zero, $zero, .L8098F040   
/* 00B48 8098EFA8 A08300C8 */  sb      $v1, 0x00C8($a0)           ## 000000C8
.L8098EFAC:
/* 00B4C 8098EFAC 248201A0 */  addiu   $v0, $a0, 0x01A0           ## $v0 = 000001A0
/* 00B50 8098EFB0 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 00B54 8098EFB4 44818000 */  mtc1    $at, $f16                  ## $f16 = 1.00
/* 00B58 8098EFB8 C44A0000 */  lwc1    $f10, 0x0000($v0)          ## 000001A0
/* 00B5C 8098EFBC 44801000 */  mtc1    $zero, $f2                 ## $f2 = 0.00
/* 00B60 8098EFC0 46105481 */  sub.s   $f18, $f10, $f16           
/* 00B64 8098EFC4 E4520000 */  swc1    $f18, 0x0000($v0)          ## 000001A0
/* 00B68 8098EFC8 C4400000 */  lwc1    $f0, 0x0000($v0)           ## 000001A0
/* 00B6C 8098EFCC 4602003E */  c.le.s  $f0, $f2                   
/* 00B70 8098EFD0 00000000 */  nop
/* 00B74 8098EFD4 45000007 */  bc1f    .L8098EFF4                 
/* 00B78 8098EFD8 240A0007 */  addiu   $t2, $zero, 0x0007         ## $t2 = 00000007
/* 00B7C 8098EFDC AC8A0198 */  sw      $t2, 0x0198($a0)           ## 00000198
/* 00B80 8098EFE0 AC80019C */  sw      $zero, 0x019C($a0)         ## 0000019C
/* 00B84 8098EFE4 E4420000 */  swc1    $f2, 0x0000($v0)           ## 000001A0
/* 00B88 8098EFE8 AC8001A4 */  sw      $zero, 0x01A4($a0)         ## 000001A4
/* 00B8C 8098EFEC 10000014 */  beq     $zero, $zero, .L8098F040   
/* 00B90 8098EFF0 A08000C8 */  sb      $zero, 0x00C8($a0)         ## 000000C8
.L8098EFF4:
/* 00B94 8098EFF4 3C058016 */  lui     $a1, 0x8016                ## $a1 = 80160000
/* 00B98 8098EFF8 24A5FA90 */  addiu   $a1, $a1, 0xFA90           ## $a1 = 8015FA90
/* 00B9C 8098EFFC 8CAB0000 */  lw      $t3, 0x0000($a1)           ## 8015FA90
/* 00BA0 8098F000 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 00BA4 8098F004 44811000 */  mtc1    $at, $f2                   ## $f2 = 10.00
/* 00BA8 8098F008 856C145E */  lh      $t4, 0x145E($t3)           ## 0000145E
/* 00BAC 8098F00C 3C01437F */  lui     $at, 0x437F                ## $at = 437F0000
/* 00BB0 8098F010 44818000 */  mtc1    $at, $f16                  ## $f16 = 255.00
/* 00BB4 8098F014 448C2000 */  mtc1    $t4, $f4                   ## $f4 = 0.00
/* 00BB8 8098F018 00000000 */  nop
/* 00BBC 8098F01C 468021A0 */  cvt.s.w $f6, $f4                   
/* 00BC0 8098F020 46023200 */  add.s   $f8, $f6, $f2              
/* 00BC4 8098F024 46080283 */  div.s   $f10, $f0, $f8             
/* 00BC8 8098F028 46105482 */  mul.s   $f18, $f10, $f16           
/* 00BCC 8098F02C 4600910D */  trunc.w.s $f4, $f18                  
/* 00BD0 8098F030 44022000 */  mfc1    $v0, $f4                   
/* 00BD4 8098F034 00000000 */  nop
/* 00BD8 8098F038 AC8201A4 */  sw      $v0, 0x01A4($a0)           ## 000001A4
/* 00BDC 8098F03C A08200C8 */  sb      $v0, 0x00C8($a0)           ## 000000C8
.L8098F040:
/* 00BE0 8098F040 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 00BE4 8098F044 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 00BE8 8098F048 03E00008 */  jr      $ra                        
/* 00BEC 8098F04C 00000000 */  nop


