glabel func_80883144
/* 00B84 80883144 27BDFFC8 */  addiu   $sp, $sp, 0xFFC8           ## $sp = FFFFFFC8
/* 00B88 80883148 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 00B8C 8088314C AFB00020 */  sw      $s0, 0x0020($sp)           
/* 00B90 80883150 AFA5003C */  sw      $a1, 0x003C($sp)           
/* 00B94 80883154 8482016A */  lh      $v0, 0x016A($a0)           ## 0000016A
/* 00B98 80883158 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00B9C 8088315C 3C014348 */  lui     $at, 0x4348                ## $at = 43480000
/* 00BA0 80883160 10400003 */  beq     $v0, $zero, .L80883170     
/* 00BA4 80883164 244EFFFF */  addiu   $t6, $v0, 0xFFFF           ## $t6 = FFFFFFFF
/* 00BA8 80883168 A48E016A */  sh      $t6, 0x016A($a0)           ## 0000016A
/* 00BAC 8088316C 8482016A */  lh      $v0, 0x016A($a0)           ## 0000016A
.L80883170:
/* 00BB0 80883170 04410004 */  bgez    $v0, .L80883184            
/* 00BB4 80883174 304F0003 */  andi    $t7, $v0, 0x0003           ## $t7 = 00000000
/* 00BB8 80883178 11E00002 */  beq     $t7, $zero, .L80883184     
/* 00BBC 8088317C 00000000 */  nop
/* 00BC0 80883180 25EFFFFC */  addiu   $t7, $t7, 0xFFFC           ## $t7 = FFFFFFFC
.L80883184:
/* 00BC4 80883184 15E00028 */  bne     $t7, $zero, .L80883228     
/* 00BC8 80883188 00000000 */  nop
/* 00BCC 8088318C 44816000 */  mtc1    $at, $f12                  ## $f12 = 200.00
/* 00BD0 80883190 0C00CFC8 */  jal     Math_Rand_CenteredFloat
              
/* 00BD4 80883194 00000000 */  nop
/* 00BD8 80883198 3C014260 */  lui     $at, 0x4260                ## $at = 42600000
/* 00BDC 8088319C 44813000 */  mtc1    $at, $f6                   ## $f6 = 56.00
/* 00BE0 808831A0 C6040024 */  lwc1    $f4, 0x0024($s0)           ## 00000024
/* 00BE4 808831A4 46062201 */  sub.s   $f8, $f4, $f6              
/* 00BE8 808831A8 46080280 */  add.s   $f10, $f0, $f8             
/* 00BEC 808831AC 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 00BF0 808831B0 E7AA002C */  swc1    $f10, 0x002C($sp)          
/* 00BF4 808831B4 3C0142A0 */  lui     $at, 0x42A0                ## $at = 42A00000
/* 00BF8 808831B8 44818000 */  mtc1    $at, $f16                  ## $f16 = 80.00
/* 00BFC 808831BC C6040028 */  lwc1    $f4, 0x0028($s0)           ## 00000028
/* 00C00 808831C0 3C014348 */  lui     $at, 0x4348                ## $at = 43480000
/* 00C04 808831C4 46100482 */  mul.s   $f18, $f0, $f16            
/* 00C08 808831C8 44816000 */  mtc1    $at, $f12                  ## $f12 = 200.00
/* 00C0C 808831CC 46049180 */  add.s   $f6, $f18, $f4             
/* 00C10 808831D0 0C00CFC8 */  jal     Math_Rand_CenteredFloat
              
/* 00C14 808831D4 E7A60030 */  swc1    $f6, 0x0030($sp)           
/* 00C18 808831D8 3C014260 */  lui     $at, 0x4260                ## $at = 42600000
/* 00C1C 808831DC 44815000 */  mtc1    $at, $f10                  ## $f10 = 56.00
/* 00C20 808831E0 C608002C */  lwc1    $f8, 0x002C($s0)           ## 0000002C
/* 00C24 808831E4 3C068088 */  lui     $a2, %hi(D_8088361C)       ## $a2 = 80880000
/* 00C28 808831E8 24C6361C */  addiu   $a2, $a2, %lo(D_8088361C)  ## $a2 = 8088361C
/* 00C2C 808831EC 460A4400 */  add.s   $f16, $f8, $f10            
/* 00C30 808831F0 24180096 */  addiu   $t8, $zero, 0x0096         ## $t8 = 00000096
/* 00C34 808831F4 24190046 */  addiu   $t9, $zero, 0x0046         ## $t9 = 00000046
/* 00C38 808831F8 AFB90014 */  sw      $t9, 0x0014($sp)           
/* 00C3C 808831FC 46100480 */  add.s   $f18, $f0, $f16            
/* 00C40 80883200 AFB80010 */  sw      $t8, 0x0010($sp)           
/* 00C44 80883204 00C03825 */  or      $a3, $a2, $zero            ## $a3 = 8088361C
/* 00C48 80883208 8FA4003C */  lw      $a0, 0x003C($sp)           
/* 00C4C 8088320C E7B20034 */  swc1    $f18, 0x0034($sp)          
/* 00C50 80883210 0C00A3A1 */  jal     func_80028E84              
/* 00C54 80883214 27A5002C */  addiu   $a1, $sp, 0x002C           ## $a1 = FFFFFFF4
/* 00C58 80883218 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00C5C 8088321C 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 00C60 80883220 2405180E */  addiu   $a1, $zero, 0x180E         ## $a1 = 0000180E
/* 00C64 80883224 8602016A */  lh      $v0, 0x016A($s0)           ## 0000016A
.L80883228:
/* 00C68 80883228 14400005 */  bne     $v0, $zero, .L80883240     
/* 00C6C 8088322C 24080014 */  addiu   $t0, $zero, 0x0014         ## $t0 = 00000014
/* 00C70 80883230 3C098088 */  lui     $t1, %hi(func_80883254)    ## $t1 = 80880000
/* 00C74 80883234 25293254 */  addiu   $t1, $t1, %lo(func_80883254) ## $t1 = 80883254
/* 00C78 80883238 A608016A */  sh      $t0, 0x016A($s0)           ## 0000016A
/* 00C7C 8088323C AE090164 */  sw      $t1, 0x0164($s0)           ## 00000164
.L80883240:
/* 00C80 80883240 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 00C84 80883244 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 00C88 80883248 27BD0038 */  addiu   $sp, $sp, 0x0038           ## $sp = 00000000
/* 00C8C 8088324C 03E00008 */  jr      $ra                        
/* 00C90 80883250 00000000 */  nop


