glabel func_809968D4
/* 00634 809968D4 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 00638 809968D8 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 0063C 809968DC AFA40028 */  sw      $a0, 0x0028($sp)           
/* 00640 809968E0 8CAE1C44 */  lw      $t6, 0x1C44($a1)           ## 00001C44
/* 00644 809968E4 AFA5002C */  sw      $a1, 0x002C($sp)           
/* 00648 809968E8 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 0064C 809968EC 0C023A62 */  jal     func_8008E988              
/* 00650 809968F0 AFAE0024 */  sw      $t6, 0x0024($sp)           
/* 00654 809968F4 14400052 */  bne     $v0, $zero, .L80996A40     
/* 00658 809968F8 8FA50028 */  lw      $a1, 0x0028($sp)           
/* 0065C 809968FC 90A2016C */  lbu     $v0, 0x016C($a1)           ## 0000016C
/* 00660 80996900 3C18809A */  lui     $t8, %hi(D_80998134)       ## $t8 = 809A0000
/* 00664 80996904 27188134 */  addiu   $t8, $t8, %lo(D_80998134)  ## $t8 = 80998134
/* 00668 80996908 00027880 */  sll     $t7, $v0,  2               
/* 0066C 8099690C 01E27823 */  subu    $t7, $t7, $v0              
/* 00670 80996910 000F7880 */  sll     $t7, $t7,  2               
/* 00674 80996914 24010003 */  addiu   $at, $zero, 0x0003         ## $at = 00000003
/* 00678 80996918 10410004 */  beq     $v0, $at, .L8099692C       
/* 0067C 8099691C 01F81821 */  addu    $v1, $t7, $t8              
/* 00680 80996920 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 00684 80996924 10000005 */  beq     $zero, $zero, .L8099693C   
/* 00688 80996928 9079000A */  lbu     $t9, 0x000A($v1)           ## 0000000A
.L8099692C:
/* 0068C 8099692C 3C0142A0 */  lui     $at, 0x42A0                ## $at = 42A00000
/* 00690 80996930 44810000 */  mtc1    $at, $f0                   ## $f0 = 80.00
/* 00694 80996934 00000000 */  nop
/* 00698 80996938 9079000A */  lbu     $t9, 0x000A($v1)           ## 0000000A
.L8099693C:
/* 0069C 8099693C 44060000 */  mfc1    $a2, $f0                   
/* 006A0 80996940 8FA4002C */  lw      $a0, 0x002C($sp)           
/* 006A4 80996944 44992000 */  mtc1    $t9, $f4                   ## $f4 = 0.00
/* 006A8 80996948 3C014F80 */  lui     $at, 0x4F80                ## $at = 4F800000
/* 006AC 8099694C 07210004 */  bgez    $t9, .L80996960            
/* 006B0 80996950 46802120 */  cvt.s.w $f4, $f4                   
/* 006B4 80996954 44813000 */  mtc1    $at, $f6                   ## $f6 = 4294967296.00
/* 006B8 80996958 00000000 */  nop
/* 006BC 8099695C 46062100 */  add.s   $f4, $f4, $f6              
.L80996960:
/* 006C0 80996960 9068000B */  lbu     $t0, 0x000B($v1)           ## 0000000B
/* 006C4 80996964 44072000 */  mfc1    $a3, $f4                   
/* 006C8 80996968 3C014F80 */  lui     $at, 0x4F80                ## $at = 4F800000
/* 006CC 8099696C 44884000 */  mtc1    $t0, $f8                   ## $f8 = 0.00
/* 006D0 80996970 05010004 */  bgez    $t0, .L80996984            
/* 006D4 80996974 468042A0 */  cvt.s.w $f10, $f8                  
/* 006D8 80996978 44818000 */  mtc1    $at, $f16                  ## $f16 = 4294967296.00
/* 006DC 8099697C 00000000 */  nop
/* 006E0 80996980 46105280 */  add.s   $f10, $f10, $f16           
.L80996984:
/* 006E4 80996984 0C265A10 */  jal     func_80996840              
/* 006E8 80996988 E7AA0010 */  swc1    $f10, 0x0010($sp)          
/* 006EC 8099698C 3C014248 */  lui     $at, 0x4248                ## $at = 42480000
/* 006F0 80996990 44819000 */  mtc1    $at, $f18                  ## $f18 = 50.00
/* 006F4 80996994 46000086 */  mov.s   $f2, $f0                   
/* 006F8 80996998 46000005 */  abs.s   $f0, $f0                   
/* 006FC 8099699C 44806000 */  mtc1    $zero, $f12                ## $f12 = 0.00
/* 00700 809969A0 4612003C */  c.lt.s  $f0, $f18                  
/* 00704 809969A4 8FA90024 */  lw      $t1, 0x0024($sp)           
/* 00708 809969A8 8FAB0028 */  lw      $t3, 0x0028($sp)           
/* 0070C 809969AC 45020025 */  bc1fl   .L80996A44                 
/* 00710 809969B0 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 00714 809969B4 852A00B6 */  lh      $t2, 0x00B6($t1)           ## 000000B6
/* 00718 809969B8 856C00B6 */  lh      $t4, 0x00B6($t3)           ## 000000B6
/* 0071C 809969BC 4602603C */  c.lt.s  $f12, $f2                  
/* 00720 809969C0 340D8000 */  ori     $t5, $zero, 0x8000         ## $t5 = 00008000
/* 00724 809969C4 014C1023 */  subu    $v0, $t2, $t4              
/* 00728 809969C8 00021400 */  sll     $v0, $v0, 16               
/* 0072C 809969CC 45000004 */  bc1f    .L809969E0                 
/* 00730 809969D0 00021403 */  sra     $v0, $v0, 16               
/* 00734 809969D4 01A21023 */  subu    $v0, $t5, $v0              
/* 00738 809969D8 00021400 */  sll     $v0, $v0, 16               
/* 0073C 809969DC 00021403 */  sra     $v0, $v0, 16               
.L809969E0:
/* 00740 809969E0 04400003 */  bltz    $v0, .L809969F0            
/* 00744 809969E4 00021823 */  subu    $v1, $zero, $v0            
/* 00748 809969E8 10000001 */  beq     $zero, $zero, .L809969F0   
/* 0074C 809969EC 00401825 */  or      $v1, $v0, $zero            ## $v1 = 00000000
.L809969F0:
/* 00750 809969F0 28613000 */  slti    $at, $v1, 0x3000           
/* 00754 809969F4 50200013 */  beql    $at, $zero, .L80996A44     
/* 00758 809969F8 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
/* 0075C 809969FC 4602603E */  c.le.s  $f12, $f2                  
/* 00760 80996A00 3C01BF80 */  lui     $at, 0xBF80                ## $at = BF800000
/* 00764 80996A04 45020009 */  bc1fl   .L80996A2C                 
/* 00768 80996A08 44810000 */  mtc1    $at, $f0                   ## $f0 = -1.00
/* 0076C 80996A0C 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 00770 80996A10 44810000 */  mtc1    $at, $f0                   ## $f0 = 1.00
/* 00774 80996A14 00000000 */  nop
/* 00778 80996A18 4600010D */  trunc.w.s $f4, $f0                   
/* 0077C 80996A1C 44022000 */  mfc1    $v0, $f4                   
/* 00780 80996A20 10000009 */  beq     $zero, $zero, .L80996A48   
/* 00784 80996A24 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 00788 80996A28 44810000 */  mtc1    $at, $f0                   ## $f0 = 1.00
.L80996A2C:
/* 0078C 80996A2C 00000000 */  nop
/* 00790 80996A30 4600010D */  trunc.w.s $f4, $f0                   
/* 00794 80996A34 44022000 */  mfc1    $v0, $f4                   
/* 00798 80996A38 10000003 */  beq     $zero, $zero, .L80996A48   
/* 0079C 80996A3C 8FBF001C */  lw      $ra, 0x001C($sp)           
.L80996A40:
/* 007A0 80996A40 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
.L80996A44:
/* 007A4 80996A44 8FBF001C */  lw      $ra, 0x001C($sp)           
.L80996A48:
/* 007A8 80996A48 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 007AC 80996A4C 03E00008 */  jr      $ra                        
/* 007B0 80996A50 00000000 */  nop


