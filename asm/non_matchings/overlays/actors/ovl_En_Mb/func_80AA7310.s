glabel func_80AA7310
/* 012C0 80AA7310 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 012C4 80AA7314 27BDFFB8 */  addiu   $sp, $sp, 0xFFB8           ## $sp = FFFFFFB8
/* 012C8 80AA7318 AFB00030 */  sw      $s0, 0x0030($sp)           
/* 012CC 80AA731C AFA5004C */  sw      $a1, 0x004C($sp)           
/* 012D0 80AA7320 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 012D4 80AA7324 AFBF0034 */  sw      $ra, 0x0034($sp)           
/* 012D8 80AA7328 44050000 */  mfc1    $a1, $f0                   
/* 012DC 80AA732C 24840068 */  addiu   $a0, $a0, 0x0068           ## $a0 = 00000068
/* 012E0 80AA7330 3C063F00 */  lui     $a2, 0x3F00                ## $a2 = 3F000000
/* 012E4 80AA7334 3C073F80 */  lui     $a3, 0x3F80                ## $a3 = 3F800000
/* 012E8 80AA7338 0C01E0C4 */  jal     Math_SmoothScaleMaxMinF
              
/* 012EC 80AA733C E7A00010 */  swc1    $f0, 0x0010($sp)           
/* 012F0 80AA7340 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 012F4 80AA7344 44812000 */  mtc1    $at, $f4                   ## $f4 = 1.00
/* 012F8 80AA7348 C6060068 */  lwc1    $f6, 0x0068($s0)           ## 00000068
/* 012FC 80AA734C 3C014080 */  lui     $at, 0x4080                ## $at = 40800000
/* 01300 80AA7350 8FA4004C */  lw      $a0, 0x004C($sp)           
/* 01304 80AA7354 4606203C */  c.lt.s  $f4, $f6                   
/* 01308 80AA7358 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 0130C 80AA735C 26060024 */  addiu   $a2, $s0, 0x0024           ## $a2 = 00000024
/* 01310 80AA7360 3C0740A0 */  lui     $a3, 0x40A0                ## $a3 = 40A00000
/* 01314 80AA7364 4500000A */  bc1f    .L80AA7390                 
/* 01318 80AA7368 240E0003 */  addiu   $t6, $zero, 0x0003         ## $t6 = 00000003
/* 0131C 80AA736C 44814000 */  mtc1    $at, $f8                   ## $f8 = 4.00
/* 01320 80AA7370 240F0064 */  addiu   $t7, $zero, 0x0064         ## $t7 = 00000064
/* 01324 80AA7374 2418000F */  addiu   $t8, $zero, 0x000F         ## $t8 = 0000000F
/* 01328 80AA7378 AFB8001C */  sw      $t8, 0x001C($sp)           
/* 0132C 80AA737C AFAF0018 */  sw      $t7, 0x0018($sp)           
/* 01330 80AA7380 AFAE0010 */  sw      $t6, 0x0010($sp)           
/* 01334 80AA7384 AFA00020 */  sw      $zero, 0x0020($sp)         
/* 01338 80AA7388 0C00CC98 */  jal     func_80033260              
/* 0133C 80AA738C E7A80014 */  swc1    $f8, 0x0014($sp)           
.L80AA7390:
/* 01340 80AA7390 2604018C */  addiu   $a0, $s0, 0x018C           ## $a0 = 0000018C
/* 01344 80AA7394 0C02927F */  jal     SkelAnime_FrameUpdateMatrix
              
/* 01348 80AA7398 AFA4003C */  sw      $a0, 0x003C($sp)           
/* 0134C 80AA739C 50400032 */  beql    $v0, $zero, .L80AA7468     
/* 01350 80AA73A0 8FBF0034 */  lw      $ra, 0x0034($sp)           
/* 01354 80AA73A4 8619032A */  lh      $t9, 0x032A($s0)           ## 0000032A
/* 01358 80AA73A8 57200021 */  bnel    $t9, $zero, .L80AA7430     
/* 0135C 80AA73AC 860D001C */  lh      $t5, 0x001C($s0)           ## 0000001C
/* 01360 80AA73B0 8608032E */  lh      $t0, 0x032E($s0)           ## 0000032E
/* 01364 80AA73B4 3C040600 */  lui     $a0, 0x0600                ## $a0 = 06000000
/* 01368 80AA73B8 2509FFFF */  addiu   $t1, $t0, 0xFFFF           ## $t1 = FFFFFFFF
/* 0136C 80AA73BC A609032E */  sh      $t1, 0x032E($s0)           ## 0000032E
/* 01370 80AA73C0 860A032E */  lh      $t2, 0x032E($s0)           ## 0000032E
/* 01374 80AA73C4 55400028 */  bnel    $t2, $zero, .L80AA7468     
/* 01378 80AA73C8 8FBF0034 */  lw      $ra, 0x0034($sp)           
/* 0137C 80AA73CC 0C028800 */  jal     SkelAnime_GetFrameCount
              
/* 01380 80AA73D0 24842C10 */  addiu   $a0, $a0, 0x2C10           ## $a0 = 06002C10
/* 01384 80AA73D4 44825000 */  mtc1    $v0, $f10                  ## $f10 = 0.00
/* 01388 80AA73D8 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 0138C 80AA73DC 3C050600 */  lui     $a1, 0x0600                ## $a1 = 06000000
/* 01390 80AA73E0 468052A0 */  cvt.s.w $f10, $f10                 
/* 01394 80AA73E4 240B0002 */  addiu   $t3, $zero, 0x0002         ## $t3 = 00000002
/* 01398 80AA73E8 AFAB0014 */  sw      $t3, 0x0014($sp)           
/* 0139C 80AA73EC 24A52C10 */  addiu   $a1, $a1, 0x2C10           ## $a1 = 06002C10
/* 013A0 80AA73F0 8FA4003C */  lw      $a0, 0x003C($sp)           
/* 013A4 80AA73F4 3C06BF80 */  lui     $a2, 0xBF80                ## $a2 = BF800000
/* 013A8 80AA73F8 44075000 */  mfc1    $a3, $f10                  
/* 013AC 80AA73FC E7A00010 */  swc1    $f0, 0x0010($sp)           
/* 013B0 80AA7400 0C029468 */  jal     SkelAnime_ChangeAnimation
              
/* 013B4 80AA7404 E7A00018 */  swc1    $f0, 0x0018($sp)           
/* 013B8 80AA7408 44808000 */  mtc1    $zero, $f16                ## $f16 = 0.00
/* 013BC 80AA740C 240C0001 */  addiu   $t4, $zero, 0x0001         ## $t4 = 00000001
/* 013C0 80AA7410 A60C032A */  sh      $t4, 0x032A($s0)           ## 0000032A
/* 013C4 80AA7414 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 013C8 80AA7418 240538BD */  addiu   $a1, $zero, 0x38BD         ## $a1 = 000038BD
/* 013CC 80AA741C 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 013D0 80AA7420 E6100068 */  swc1    $f16, 0x0068($s0)          ## 00000068
/* 013D4 80AA7424 10000010 */  beq     $zero, $zero, .L80AA7468   
/* 013D8 80AA7428 8FBF0034 */  lw      $ra, 0x0034($sp)           
/* 013DC 80AA742C 860D001C */  lh      $t5, 0x001C($s0)           ## 0000001C
.L80AA7430:
/* 013E0 80AA7430 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 013E4 80AA7434 05A10009 */  bgez    $t5, .L80AA745C            
/* 013E8 80AA7438 00000000 */  nop
/* 013EC 80AA743C 0C2A9A5D */  jal     func_80AA6974              
/* 013F0 80AA7440 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 013F4 80AA7444 240E0050 */  addiu   $t6, $zero, 0x0050         ## $t6 = 00000050
/* 013F8 80AA7448 A60E032E */  sh      $t6, 0x032E($s0)           ## 0000032E
/* 013FC 80AA744C 8602032E */  lh      $v0, 0x032E($s0)           ## 0000032E
/* 01400 80AA7450 A602032C */  sh      $v0, 0x032C($s0)           ## 0000032C
/* 01404 80AA7454 10000003 */  beq     $zero, $zero, .L80AA7464   
/* 01408 80AA7458 A602032A */  sh      $v0, 0x032A($s0)           ## 0000032A
.L80AA745C:
/* 0140C 80AA745C 0C2A9A3F */  jal     func_80AA68FC              
/* 01410 80AA7460 8FA5004C */  lw      $a1, 0x004C($sp)           
.L80AA7464:
/* 01414 80AA7464 8FBF0034 */  lw      $ra, 0x0034($sp)           
.L80AA7468:
/* 01418 80AA7468 8FB00030 */  lw      $s0, 0x0030($sp)           
/* 0141C 80AA746C 27BD0048 */  addiu   $sp, $sp, 0x0048           ## $sp = 00000000
/* 01420 80AA7470 03E00008 */  jr      $ra                        
/* 01424 80AA7474 00000000 */  nop


