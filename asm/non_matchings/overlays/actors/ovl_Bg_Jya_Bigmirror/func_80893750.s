glabel func_80893750
/* 00070 80893750 27BDFFA0 */  addiu   $sp, $sp, 0xFFA0           ## $sp = FFFFFFA0
/* 00074 80893754 AFBF005C */  sw      $ra, 0x005C($sp)           
/* 00078 80893758 AFBE0058 */  sw      $s8, 0x0058($sp)           
/* 0007C 8089375C AFB70054 */  sw      $s7, 0x0054($sp)           
/* 00080 80893760 AFB60050 */  sw      $s6, 0x0050($sp)           
/* 00084 80893764 AFB5004C */  sw      $s5, 0x004C($sp)           
/* 00088 80893768 AFB40048 */  sw      $s4, 0x0048($sp)           
/* 0008C 8089376C AFB30044 */  sw      $s3, 0x0044($sp)           
/* 00090 80893770 AFB20040 */  sw      $s2, 0x0040($sp)           
/* 00094 80893774 AFB1003C */  sw      $s1, 0x003C($sp)           
/* 00098 80893778 AFB00038 */  sw      $s0, 0x0038($sp)           
/* 0009C 8089377C 908E015C */  lbu     $t6, 0x015C($a0)           ## 0000015C
/* 000A0 80893780 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 000A4 80893784 00A0A825 */  or      $s5, $a1, $zero            ## $s5 = 00000000
/* 000A8 80893788 31CF0030 */  andi    $t7, $t6, 0x0030           ## $t7 = 00000000
/* 000AC 8089378C 11E0004E */  beq     $t7, $zero, .L808938C8     
/* 000B0 80893790 00009025 */  or      $s2, $zero, $zero          ## $s2 = 00000000
/* 000B4 80893794 3C128089 */  lui     $s2, %hi(D_80893EF4)       ## $s2 = 80890000
/* 000B8 80893798 3C178089 */  lui     $s7, %hi(D_80893F84)       ## $s7 = 80890000
/* 000BC 8089379C 3C168089 */  lui     $s6, %hi(D_80893F60)       ## $s6 = 80890000
/* 000C0 808937A0 26D63F60 */  addiu   $s6, $s6, %lo(D_80893F60)  ## $s6 = 80893F60
/* 000C4 808937A4 26F73F84 */  addiu   $s7, $s7, %lo(D_80893F84)  ## $s7 = 80893F84
/* 000C8 808937A8 26523EF4 */  addiu   $s2, $s2, %lo(D_80893EF4)  ## $s2 = 80893EF4
/* 000CC 808937AC 0000A025 */  or      $s4, $zero, $zero          ## $s4 = 00000000
/* 000D0 808937B0 00809825 */  or      $s3, $a0, $zero            ## $s3 = 00000000
/* 000D4 808937B4 2491014C */  addiu   $s1, $a0, 0x014C           ## $s1 = 0000014C
/* 000D8 808937B8 241E0002 */  addiu   $s8, $zero, 0x0002         ## $s8 = 00000002
.L808937BC:
/* 000DC 808937BC 8E78014C */  lw      $t8, 0x014C($s3)           ## 0000014C
/* 000E0 808937C0 53000020 */  beql    $t8, $zero, .L80893844     
/* 000E4 808937C4 C6440000 */  lwc1    $f4, 0x0000($s2)           ## 80893EF4
/* 000E8 808937C8 8E390000 */  lw      $t9, 0x0000($s1)           ## 0000014C
/* 000EC 808937CC 3C098089 */  lui     $t1, %hi(D_80893F1C)       ## $t1 = 80890000
/* 000F0 808937D0 25293F1C */  addiu   $t1, $t1, %lo(D_80893F1C)  ## $t1 = 80893F1C
/* 000F4 808937D4 872800B6 */  lh      $t0, 0x00B6($t9)           ## 000000B6
/* 000F8 808937D8 02891021 */  addu    $v0, $s4, $t1              
/* 000FC 808937DC 02C02025 */  or      $a0, $s6, $zero            ## $a0 = 80893F60
/* 00100 808937E0 A6280004 */  sh      $t0, 0x0004($s1)           ## 00000150
/* 00104 808937E4 862B0004 */  lh      $t3, 0x0004($s1)           ## 00000150
/* 00108 808937E8 864A000E */  lh      $t2, 0x000E($s2)           ## 80893F02
/* 0010C 808937EC 02E02825 */  or      $a1, $s7, $zero            ## $a1 = 80893F84
/* 00110 808937F0 554B0007 */  bnel    $t2, $t3, .L80893810       
/* 00114 808937F4 90580000 */  lbu     $t8, 0x0000($v0)           ## 00000000
/* 00118 808937F8 920C015C */  lbu     $t4, 0x015C($s0)           ## 0000015C
/* 0011C 808937FC 904D0000 */  lbu     $t5, 0x0000($v0)           ## 00000000
/* 00120 80893800 018D7025 */  or      $t6, $t4, $t5              ## $t6 = 00000000
/* 00124 80893804 10000006 */  beq     $zero, $zero, .L80893820   
/* 00128 80893808 A20E015C */  sb      $t6, 0x015C($s0)           ## 0000015C
/* 0012C 8089380C 90580000 */  lbu     $t8, 0x0000($v0)           ## 00000000
.L80893810:
/* 00130 80893810 920F015C */  lbu     $t7, 0x015C($s0)           ## 0000015C
/* 00134 80893814 0300C827 */  nor     $t9, $t8, $zero            
/* 00138 80893818 01F94024 */  and     $t0, $t7, $t9              
/* 0013C 8089381C A208015C */  sb      $t0, 0x015C($s0)           ## 0000015C
.L80893820:
/* 00140 80893820 8E290000 */  lw      $t1, 0x0000($s1)           ## 0000014C
/* 00144 80893824 8D2A0130 */  lw      $t2, 0x0130($t1)           ## 8089404C
/* 00148 80893828 55400021 */  bnel    $t2, $zero, .L808938B0     
/* 0014C 8089382C 26940001 */  addiu   $s4, $s4, 0x0001           ## $s4 = 00000001
/* 00150 80893830 0C00084C */  jal     osSyncPrintf
              
/* 00154 80893834 240600CB */  addiu   $a2, $zero, 0x00CB         ## $a2 = 000000CB
/* 00158 80893838 1000001D */  beq     $zero, $zero, .L808938B0   
/* 0015C 8089383C 26940001 */  addiu   $s4, $s4, 0x0001           ## $s4 = 00000002
/* 00160 80893840 C6440000 */  lwc1    $f4, 0x0000($s2)           ## 80893EF4
.L80893844:
/* 00164 80893844 26A41C24 */  addiu   $a0, $s5, 0x1C24           ## $a0 = 00001C24
/* 00168 80893848 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0016C 8089384C E7A40010 */  swc1    $f4, 0x0010($sp)           
/* 00170 80893850 C6460004 */  lwc1    $f6, 0x0004($s2)           ## 80893EF8
/* 00174 80893854 02A03025 */  or      $a2, $s5, $zero            ## $a2 = 00000000
/* 00178 80893858 240700FC */  addiu   $a3, $zero, 0x00FC         ## $a3 = 000000FC
/* 0017C 8089385C E7A60014 */  swc1    $f6, 0x0014($sp)           
/* 00180 80893860 C6480008 */  lwc1    $f8, 0x0008($s2)           ## 80893EFC
/* 00184 80893864 AFA0001C */  sw      $zero, 0x001C($sp)         
/* 00188 80893868 E7A80018 */  swc1    $f8, 0x0018($sp)           
/* 0018C 8089386C 862B0004 */  lh      $t3, 0x0004($s1)           ## 00000150
/* 00190 80893870 AFA00024 */  sw      $zero, 0x0024($sp)         
/* 00194 80893874 AFAB0020 */  sw      $t3, 0x0020($sp)           
/* 00198 80893878 864C000C */  lh      $t4, 0x000C($s2)           ## 80893F00
/* 0019C 8089387C 0C00C916 */  jal     Actor_SpawnAttached
              
/* 001A0 80893880 AFAC0028 */  sw      $t4, 0x0028($sp)           
/* 001A4 80893884 AE220000 */  sw      $v0, 0x0000($s1)           ## 0000014C
/* 001A8 80893888 AE00011C */  sw      $zero, 0x011C($s0)         ## 0000011C
/* 001AC 8089388C 8E2D0000 */  lw      $t5, 0x0000($s1)           ## 0000014C
/* 001B0 80893890 3C048089 */  lui     $a0, %hi(D_80893F9C)       ## $a0 = 80890000
/* 001B4 80893894 24843F9C */  addiu   $a0, $a0, %lo(D_80893F9C)  ## $a0 = 80893F9C
/* 001B8 80893898 15A00004 */  bne     $t5, $zero, .L808938AC     
/* 001BC 8089389C 3C058089 */  lui     $a1, %hi(D_80893FBC)       ## $a1 = 80890000
/* 001C0 808938A0 24A53FBC */  addiu   $a1, $a1, %lo(D_80893FBC)  ## $a1 = 80893FBC
/* 001C4 808938A4 0C00084C */  jal     osSyncPrintf
              
/* 001C8 808938A8 240600DD */  addiu   $a2, $zero, 0x00DD         ## $a2 = 000000DD
.L808938AC:
/* 001CC 808938AC 26940001 */  addiu   $s4, $s4, 0x0001           ## $s4 = 00000003
.L808938B0:
/* 001D0 808938B0 26520014 */  addiu   $s2, $s2, 0x0014           ## $s2 = 80893F08
/* 001D4 808938B4 26730008 */  addiu   $s3, $s3, 0x0008           ## $s3 = 00000008
/* 001D8 808938B8 169EFFC0 */  bne     $s4, $s8, .L808937BC       
/* 001DC 808938BC 26310008 */  addiu   $s1, $s1, 0x0008           ## $s1 = 00000154
/* 001E0 808938C0 10000017 */  beq     $zero, $zero, .L80893920   
/* 001E4 808938C4 8FBF005C */  lw      $ra, 0x005C($sp)           
.L808938C8:
/* 001E8 808938C8 02009825 */  or      $s3, $s0, $zero            ## $s3 = 00000000
/* 001EC 808938CC 24140010 */  addiu   $s4, $zero, 0x0010         ## $s4 = 00000010
.L808938D0:
/* 001F0 808938D0 8E6E014C */  lw      $t6, 0x014C($s3)           ## 0000014C
/* 001F4 808938D4 2671014C */  addiu   $s1, $s3, 0x014C           ## $s1 = 0000014C
/* 001F8 808938D8 51C0000E */  beql    $t6, $zero, .L80893914     
/* 001FC 808938DC 26520008 */  addiu   $s2, $s2, 0x0008           ## $s2 = 80893F10
/* 00200 808938E0 8E240000 */  lw      $a0, 0x0000($s1)           ## 0000014C
/* 00204 808938E4 8C85011C */  lw      $a1, 0x011C($a0)           ## 0000011C
/* 00208 808938E8 10A00006 */  beq     $a1, $zero, .L80893904     
/* 0020C 808938EC 00000000 */  nop
/* 00210 808938F0 0C00B55C */  jal     Actor_Kill
              
/* 00214 808938F4 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 00218 808938F8 8E380000 */  lw      $t8, 0x0000($s1)           ## 0000014C
/* 0021C 808938FC AF00011C */  sw      $zero, 0x011C($t8)         ## 0000011C
/* 00220 80893900 8E240000 */  lw      $a0, 0x0000($s1)           ## 0000014C
.L80893904:
/* 00224 80893904 0C00B55C */  jal     Actor_Kill
              
/* 00228 80893908 00000000 */  nop
/* 0022C 8089390C AE200000 */  sw      $zero, 0x0000($s1)         ## 0000014C
/* 00230 80893910 26520008 */  addiu   $s2, $s2, 0x0008           ## $s2 = 80893F18
.L80893914:
/* 00234 80893914 1654FFEE */  bne     $s2, $s4, .L808938D0       
/* 00238 80893918 26730008 */  addiu   $s3, $s3, 0x0008           ## $s3 = 00000008
/* 0023C 8089391C 8FBF005C */  lw      $ra, 0x005C($sp)           
.L80893920:
/* 00240 80893920 8FB00038 */  lw      $s0, 0x0038($sp)           
/* 00244 80893924 8FB1003C */  lw      $s1, 0x003C($sp)           
/* 00248 80893928 8FB20040 */  lw      $s2, 0x0040($sp)           
/* 0024C 8089392C 8FB30044 */  lw      $s3, 0x0044($sp)           
/* 00250 80893930 8FB40048 */  lw      $s4, 0x0048($sp)           
/* 00254 80893934 8FB5004C */  lw      $s5, 0x004C($sp)           
/* 00258 80893938 8FB60050 */  lw      $s6, 0x0050($sp)           
/* 0025C 8089393C 8FB70054 */  lw      $s7, 0x0054($sp)           
/* 00260 80893940 8FBE0058 */  lw      $s8, 0x0058($sp)           
/* 00264 80893944 03E00008 */  jr      $ra                        
/* 00268 80893948 27BD0060 */  addiu   $sp, $sp, 0x0060           ## $sp = 00000000


