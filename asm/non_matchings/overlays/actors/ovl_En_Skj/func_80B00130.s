glabel func_80B00130
/* 01E80 80B00130 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 01E84 80B00134 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 01E88 80B00138 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 01E8C 80B0013C AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 01E90 80B00140 AFA5002C */  sw      $a1, 0x002C($sp)           
/* 01E94 80B00144 C60402F0 */  lwc1    $f4, 0x02F0($s0)           ## 000002F0
/* 01E98 80B00148 8605008A */  lh      $a1, 0x008A($s0)           ## 0000008A
/* 01E9C 80B0014C AFA00010 */  sw      $zero, 0x0010($sp)         
/* 01EA0 80B00150 4600218D */  trunc.w.s $f6, $f4                   
/* 01EA4 80B00154 248400B6 */  addiu   $a0, $a0, 0x00B6           ## $a0 = 000000B6
/* 01EA8 80B00158 2406000A */  addiu   $a2, $zero, 0x000A         ## $a2 = 0000000A
/* 01EAC 80B0015C 44073000 */  mfc1    $a3, $f6                   
/* 01EB0 80B00160 00000000 */  nop
/* 01EB4 80B00164 00073C00 */  sll     $a3, $a3, 16               
/* 01EB8 80B00168 0C01E1A7 */  jal     Math_SmoothScaleMaxMinS
              
/* 01EBC 80B0016C 00073C03 */  sra     $a3, $a3, 16               
/* 01EC0 80B00170 260402F0 */  addiu   $a0, $s0, 0x02F0           ## $a0 = 000002F0
/* 01EC4 80B00174 3C0544FA */  lui     $a1, 0x44FA                ## $a1 = 44FA0000
/* 01EC8 80B00178 3C063F80 */  lui     $a2, 0x3F80                ## $a2 = 3F800000
/* 01ECC 80B0017C 0C01E107 */  jal     Math_SmoothScaleMaxF
              
/* 01ED0 80B00180 3C0742C8 */  lui     $a3, 0x42C8                ## $a3 = 42C80000
/* 01ED4 80B00184 3C0142F0 */  lui     $at, 0x42F0                ## $at = 42F00000
/* 01ED8 80B00188 44815000 */  mtc1    $at, $f10                  ## $f10 = 120.00
/* 01EDC 80B0018C C6080090 */  lwc1    $f8, 0x0090($s0)           ## 00000090
/* 01EE0 80B00190 860F00B6 */  lh      $t7, 0x00B6($s0)           ## 000000B6
/* 01EE4 80B00194 460A403C */  c.lt.s  $f8, $f10                  
/* 01EE8 80B00198 A60F0032 */  sh      $t7, 0x0032($s0)           ## 00000032
/* 01EEC 80B0019C 45020007 */  bc1fl   .L80B001BC                 
/* 01EF0 80B001A0 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 01EF4 80B001A4 44808000 */  mtc1    $zero, $f16                ## $f16 = 0.00
/* 01EF8 80B001A8 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01EFC 80B001AC E6100068 */  swc1    $f16, 0x0068($s0)          ## 00000068
/* 01F00 80B001B0 0C2C0073 */  jal     func_80B001CC              
/* 01F04 80B001B4 8FA5002C */  lw      $a1, 0x002C($sp)           
/* 01F08 80B001B8 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80B001BC:
/* 01F0C 80B001BC 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 01F10 80B001C0 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 01F14 80B001C4 03E00008 */  jr      $ra                        
/* 01F18 80B001C8 00000000 */  nop


