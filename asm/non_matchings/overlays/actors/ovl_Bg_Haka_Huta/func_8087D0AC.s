glabel func_8087D0AC
/* 000EC 8087D0AC 27BDFF58 */  addiu   $sp, $sp, 0xFF58           ## $sp = FFFFFF58
/* 000F0 8087D0B0 AFBF0074 */  sw      $ra, 0x0074($sp)           
/* 000F4 8087D0B4 AFBE0070 */  sw      $s8, 0x0070($sp)           
/* 000F8 8087D0B8 AFB7006C */  sw      $s7, 0x006C($sp)           
/* 000FC 8087D0BC AFB60068 */  sw      $s6, 0x0068($sp)           
/* 00100 8087D0C0 AFB50064 */  sw      $s5, 0x0064($sp)           
/* 00104 8087D0C4 AFB40060 */  sw      $s4, 0x0060($sp)           
/* 00108 8087D0C8 AFB3005C */  sw      $s3, 0x005C($sp)           
/* 0010C 8087D0CC AFB20058 */  sw      $s2, 0x0058($sp)           
/* 00110 8087D0D0 AFB10054 */  sw      $s1, 0x0054($sp)           
/* 00114 8087D0D4 AFB00050 */  sw      $s0, 0x0050($sp)           
/* 00118 8087D0D8 F7BC0048 */  sdc1    $f28, 0x0048($sp)          
/* 0011C 8087D0DC F7BA0040 */  sdc1    $f26, 0x0040($sp)          
/* 00120 8087D0E0 F7B80038 */  sdc1    $f24, 0x0038($sp)          
/* 00124 8087D0E4 F7B60030 */  sdc1    $f22, 0x0030($sp)          
/* 00128 8087D0E8 F7B40028 */  sdc1    $f20, 0x0028($sp)          
/* 0012C 8087D0EC AFA500AC */  sw      $a1, 0x00AC($sp)           
/* 00130 8087D0F0 848E0032 */  lh      $t6, 0x0032($a0)           ## 00000032
/* 00134 8087D0F4 00808825 */  or      $s1, $a0, $zero            ## $s1 = 00000000
/* 00138 8087D0F8 3C01BF80 */  lui     $at, 0xBF80                ## $at = BF800000
/* 0013C 8087D0FC 55C00006 */  bnel    $t6, $zero, .L8087D118     
/* 00140 8087D100 4481A000 */  mtc1    $at, $f20                  ## $f20 = -1.00
/* 00144 8087D104 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 00148 8087D108 4481A000 */  mtc1    $at, $f20                  ## $f20 = 1.00
/* 0014C 8087D10C 10000004 */  beq     $zero, $zero, .L8087D120   
/* 00150 8087D110 3C01BF00 */  lui     $at, 0xBF00                ## $at = BF000000
/* 00154 8087D114 4481A000 */  mtc1    $at, $f20                  ## $f20 = -0.50
.L8087D118:
/* 00158 8087D118 00000000 */  nop
/* 0015C 8087D11C 3C01BF00 */  lui     $at, 0xBF00                ## $at = BF000000
.L8087D120:
/* 00160 8087D120 44812000 */  mtc1    $at, $f4                   ## $f4 = -0.50
/* 00164 8087D124 44801000 */  mtc1    $zero, $f2                 ## $f2 = 0.00
/* 00168 8087D128 3C014248 */  lui     $at, 0x4248                ## $at = 42480000
/* 0016C 8087D12C 46142182 */  mul.s   $f6, $f4, $f20             
/* 00170 8087D130 E7A2008C */  swc1    $f2, 0x008C($sp)           
/* 00174 8087D134 E7A20090 */  swc1    $f2, 0x0090($sp)           
/* 00178 8087D138 4481E000 */  mtc1    $at, $f28                  ## $f28 = 50.00
/* 0017C 8087D13C 3C178088 */  lui     $s7, %hi(D_8087D954)       ## $s7 = 80880000
/* 00180 8087D140 3C168088 */  lui     $s6, %hi(D_8087D950)       ## $s6 = 80880000
/* 00184 8087D144 3C158088 */  lui     $s5, %hi(D_8087D944)       ## $s5 = 80880000
/* 00188 8087D148 E7A60088 */  swc1    $f6, 0x0088($sp)           
/* 0018C 8087D14C C6280028 */  lwc1    $f8, 0x0028($s1)           ## 00000028
/* 00190 8087D150 26B5D944 */  addiu   $s5, $s5, %lo(D_8087D944)  ## $s5 = 8087D944
/* 00194 8087D154 26D6D950 */  addiu   $s6, $s6, %lo(D_8087D950)  ## $s6 = 8087D950
/* 00198 8087D158 E7A80098 */  swc1    $f8, 0x0098($sp)           
/* 0019C 8087D15C C62A002C */  lwc1    $f10, 0x002C($s1)          ## 0000002C
/* 001A0 8087D160 44814000 */  mtc1    $at, $f8                   ## $f8 = 50.00
/* 001A4 8087D164 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 001A8 8087D168 E7AA009C */  swc1    $f10, 0x009C($sp)          
/* 001AC 8087D16C C6320008 */  lwc1    $f18, 0x0008($s1)          ## 00000008
/* 001B0 8087D170 C6300024 */  lwc1    $f16, 0x0024($s1)          ## 00000024
/* 001B4 8087D174 4481D000 */  mtc1    $at, $f26                  ## $f26 = 10.00
/* 001B8 8087D178 3C0142F0 */  lui     $at, 0x42F0                ## $at = 42F00000
/* 001BC 8087D17C 46128101 */  sub.s   $f4, $f16, $f18            
/* 001C0 8087D180 4481C000 */  mtc1    $at, $f24                  ## $f24 = 120.00
/* 001C4 8087D184 26F7D954 */  addiu   $s7, $s7, %lo(D_8087D954)  ## $s7 = 8087D954
/* 001C8 8087D188 00008025 */  or      $s0, $zero, $zero          ## $s0 = 00000000
/* 001CC 8087D18C 46142182 */  mul.s   $f6, $f4, $f20             
/* 001D0 8087D190 241E0004 */  addiu   $s8, $zero, 0x0004         ## $s8 = 00000004
/* 001D4 8087D194 27B40088 */  addiu   $s4, $sp, 0x0088           ## $s4 = FFFFFFE0
/* 001D8 8087D198 27B30094 */  addiu   $s3, $sp, 0x0094           ## $s3 = FFFFFFEC
/* 001DC 8087D19C 24120002 */  addiu   $s2, $zero, 0x0002         ## $s2 = 00000002
/* 001E0 8087D1A0 46064001 */  sub.s   $f0, $f8, $f6              
/* 001E4 8087D1A4 46140582 */  mul.s   $f22, $f0, $f20            
/* 001E8 8087D1A8 00000000 */  nop
.L8087D1AC:
/* 001EC 8087D1AC 16120005 */  bne     $s0, $s2, .L8087D1C4       
/* 001F0 8087D1B0 00000000 */  nop
/* 001F4 8087D1B4 4614C402 */  mul.s   $f16, $f24, $f20           
/* 001F8 8087D1B8 C7AA009C */  lwc1    $f10, 0x009C($sp)          
/* 001FC 8087D1BC 46105480 */  add.s   $f18, $f10, $f16           
/* 00200 8087D1C0 E7B2009C */  swc1    $f18, 0x009C($sp)          
.L8087D1C4:
/* 00204 8087D1C4 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 00208 8087D1C8 00000000 */  nop
/* 0020C 8087D1CC 46160202 */  mul.s   $f8, $f0, $f22             
/* 00210 8087D1D0 C6240008 */  lwc1    $f4, 0x0008($s1)           ## 00000008
/* 00214 8087D1D4 46082181 */  sub.s   $f6, $f4, $f8              
/* 00218 8087D1D8 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 0021C 8087D1DC E7A60094 */  swc1    $f6, 0x0094($sp)           
/* 00220 8087D1E0 461A0282 */  mul.s   $f10, $f0, $f26            
/* 00224 8087D1E4 2419000A */  addiu   $t9, $zero, 0x000A         ## $t9 = 0000000A
/* 00228 8087D1E8 AFB9001C */  sw      $t9, 0x001C($sp)           
/* 0022C 8087D1EC 8FA400AC */  lw      $a0, 0x00AC($sp)           
/* 00230 8087D1F0 02602825 */  or      $a1, $s3, $zero            ## $a1 = FFFFFFEC
/* 00234 8087D1F4 02803025 */  or      $a2, $s4, $zero            ## $a2 = FFFFFFE0
/* 00238 8087D1F8 02A03825 */  or      $a3, $s5, $zero            ## $a3 = 8087D944
/* 0023C 8087D1FC 461C5400 */  add.s   $f16, $f10, $f28           
/* 00240 8087D200 AFB60010 */  sw      $s6, 0x0010($sp)           
/* 00244 8087D204 AFB70014 */  sw      $s7, 0x0014($sp)           
/* 00248 8087D208 4600848D */  trunc.w.s $f18, $f16                 
/* 0024C 8087D20C 44189000 */  mfc1    $t8, $f18                  
/* 00250 8087D210 0C00A0A7 */  jal     func_8002829C              
/* 00254 8087D214 AFB80018 */  sw      $t8, 0x0018($sp)           
/* 00258 8087D218 26100001 */  addiu   $s0, $s0, 0x0001           ## $s0 = 00000001
/* 0025C 8087D21C 161EFFE3 */  bne     $s0, $s8, .L8087D1AC       
/* 00260 8087D220 00000000 */  nop
/* 00264 8087D224 8FBF0074 */  lw      $ra, 0x0074($sp)           
/* 00268 8087D228 D7B40028 */  ldc1    $f20, 0x0028($sp)          
/* 0026C 8087D22C D7B60030 */  ldc1    $f22, 0x0030($sp)          
/* 00270 8087D230 D7B80038 */  ldc1    $f24, 0x0038($sp)          
/* 00274 8087D234 D7BA0040 */  ldc1    $f26, 0x0040($sp)          
/* 00278 8087D238 D7BC0048 */  ldc1    $f28, 0x0048($sp)          
/* 0027C 8087D23C 8FB00050 */  lw      $s0, 0x0050($sp)           
/* 00280 8087D240 8FB10054 */  lw      $s1, 0x0054($sp)           
/* 00284 8087D244 8FB20058 */  lw      $s2, 0x0058($sp)           
/* 00288 8087D248 8FB3005C */  lw      $s3, 0x005C($sp)           
/* 0028C 8087D24C 8FB40060 */  lw      $s4, 0x0060($sp)           
/* 00290 8087D250 8FB50064 */  lw      $s5, 0x0064($sp)           
/* 00294 8087D254 8FB60068 */  lw      $s6, 0x0068($sp)           
/* 00298 8087D258 8FB7006C */  lw      $s7, 0x006C($sp)           
/* 0029C 8087D25C 8FBE0070 */  lw      $s8, 0x0070($sp)           
/* 002A0 8087D260 03E00008 */  jr      $ra                        
/* 002A4 8087D264 27BD00A8 */  addiu   $sp, $sp, 0x00A8           ## $sp = 00000000


