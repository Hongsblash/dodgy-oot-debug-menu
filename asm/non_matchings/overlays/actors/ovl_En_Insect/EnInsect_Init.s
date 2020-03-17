glabel EnInsect_Init
/* 002CC 80A7C0EC 27BDFFB8 */  addiu   $sp, $sp, 0xFFB8           ## $sp = FFFFFFB8
/* 002D0 80A7C0F0 AFB3003C */  sw      $s3, 0x003C($sp)           
/* 002D4 80A7C0F4 00A09825 */  or      $s3, $a1, $zero            ## $s3 = 00000000
/* 002D8 80A7C0F8 AFBF0044 */  sw      $ra, 0x0044($sp)           
/* 002DC 80A7C0FC AFB10034 */  sw      $s1, 0x0034($sp)           
/* 002E0 80A7C100 3C0580A8 */  lui     $a1, %hi(D_80A7DF18)       ## $a1 = 80A80000
/* 002E4 80A7C104 00808825 */  or      $s1, $a0, $zero            ## $s1 = 00000000
/* 002E8 80A7C108 AFB40040 */  sw      $s4, 0x0040($sp)           
/* 002EC 80A7C10C AFB20038 */  sw      $s2, 0x0038($sp)           
/* 002F0 80A7C110 AFB00030 */  sw      $s0, 0x0030($sp)           
/* 002F4 80A7C114 0C01E037 */  jal     Actor_ProcessInitChain
              
/* 002F8 80A7C118 24A5DF18 */  addiu   $a1, $a1, %lo(D_80A7DF18)  ## $a1 = 80A7DF18
/* 002FC 80A7C11C 0C29EF88 */  jal     func_80A7BE20              
/* 00300 80A7C120 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 00304 80A7C124 8632001C */  lh      $s2, 0x001C($s1)           ## 0000001C
/* 00308 80A7C128 3C060403 */  lui     $a2, 0x0403                ## $a2 = 04030000
/* 0030C 80A7C12C 3C070403 */  lui     $a3, 0x0403                ## $a3 = 04030000
/* 00310 80A7C130 32520003 */  andi    $s2, $s2, 0x0003           ## $s2 = 00000000
/* 00314 80A7C134 00129400 */  sll     $s2, $s2, 16               
/* 00318 80A7C138 262E01F0 */  addiu   $t6, $s1, 0x01F0           ## $t6 = 000001F0
/* 0031C 80A7C13C 262F0280 */  addiu   $t7, $s1, 0x0280           ## $t7 = 00000280
/* 00320 80A7C140 24180018 */  addiu   $t8, $zero, 0x0018         ## $t8 = 00000018
/* 00324 80A7C144 00129403 */  sra     $s2, $s2, 16               
/* 00328 80A7C148 AFB80018 */  sw      $t8, 0x0018($sp)           
/* 0032C 80A7C14C AFAF0014 */  sw      $t7, 0x0014($sp)           
/* 00330 80A7C150 AFAE0010 */  sw      $t6, 0x0010($sp)           
/* 00334 80A7C154 24E741FC */  addiu   $a3, $a3, 0x41FC           ## $a3 = 040341FC
/* 00338 80A7C158 24C65590 */  addiu   $a2, $a2, 0x5590           ## $a2 = 04035590
/* 0033C 80A7C15C 02602025 */  or      $a0, $s3, $zero            ## $a0 = 00000000
/* 00340 80A7C160 0C02915F */  jal     SkelAnime_Init
              
/* 00344 80A7C164 262501AC */  addiu   $a1, $s1, 0x01AC           ## $a1 = 000001AC
/* 00348 80A7C168 2630014C */  addiu   $s0, $s1, 0x014C           ## $s0 = 0000014C
/* 0034C 80A7C16C 02002825 */  or      $a1, $s0, $zero            ## $a1 = 0000014C
/* 00350 80A7C170 0C016EFE */  jal     func_8005BBF8              
/* 00354 80A7C174 02602025 */  or      $a0, $s3, $zero            ## $a0 = 00000000
/* 00358 80A7C178 3C0780A8 */  lui     $a3, %hi(D_80A7DF00)       ## $a3 = 80A80000
/* 0035C 80A7C17C 2639016C */  addiu   $t9, $s1, 0x016C           ## $t9 = 0000016C
/* 00360 80A7C180 AFB90010 */  sw      $t9, 0x0010($sp)           
/* 00364 80A7C184 24E7DF00 */  addiu   $a3, $a3, %lo(D_80A7DF00)  ## $a3 = 80A7DF00
/* 00368 80A7C188 02602025 */  or      $a0, $s3, $zero            ## $a0 = 00000000
/* 0036C 80A7C18C 02002825 */  or      $a1, $s0, $zero            ## $a1 = 0000014C
/* 00370 80A7C190 0C017014 */  jal     func_8005C050              
/* 00374 80A7C194 02203025 */  or      $a2, $s1, $zero            ## $a2 = 00000000
/* 00378 80A7C198 96230314 */  lhu     $v1, 0x0314($s1)           ## 00000314
/* 0037C 80A7C19C 2408001E */  addiu   $t0, $zero, 0x001E         ## $t0 = 0000001E
/* 00380 80A7C1A0 A22800AE */  sb      $t0, 0x00AE($s1)           ## 000000AE
/* 00384 80A7C1A4 30690001 */  andi    $t1, $v1, 0x0001           ## $t1 = 00000000
/* 00388 80A7C1A8 11200007 */  beq     $t1, $zero, .L80A7C1C8     
/* 0038C 80A7C1AC 3C0180A8 */  lui     $at, %hi(D_80A7DFD0)       ## $at = 80A80000
/* 00390 80A7C1B0 C424DFD0 */  lwc1    $f4, %lo(D_80A7DFD0)($at)  
/* 00394 80A7C1B4 3C01C000 */  lui     $at, 0xC000                ## $at = C0000000
/* 00398 80A7C1B8 44813000 */  mtc1    $at, $f6                   ## $f6 = -2.00
/* 0039C 80A7C1BC 96230314 */  lhu     $v1, 0x0314($s1)           ## 00000314
/* 003A0 80A7C1C0 E624006C */  swc1    $f4, 0x006C($s1)           ## 0000006C
/* 003A4 80A7C1C4 E6260070 */  swc1    $f6, 0x0070($s1)           ## 00000070
.L80A7C1C8:
/* 003A8 80A7C1C8 306A0004 */  andi    $t2, $v1, 0x0004           ## $t2 = 00000000
/* 003AC 80A7C1CC 11400007 */  beq     $t2, $zero, .L80A7C1EC     
/* 003B0 80A7C1D0 240400C8 */  addiu   $a0, $zero, 0x00C8         ## $a0 = 000000C8
/* 003B4 80A7C1D4 0C01DF64 */  jal     Math_Rand_S16Offset
              
/* 003B8 80A7C1D8 24050028 */  addiu   $a1, $zero, 0x0028         ## $a1 = 00000028
/* 003BC 80A7C1DC 8E2B0004 */  lw      $t3, 0x0004($s1)           ## 00000004
/* 003C0 80A7C1E0 A622031C */  sh      $v0, 0x031C($s1)           ## 0000031C
/* 003C4 80A7C1E4 356C0010 */  ori     $t4, $t3, 0x0010           ## $t4 = 00000010
/* 003C8 80A7C1E8 AE2C0004 */  sw      $t4, 0x0004($s1)           ## 00000004
.L80A7C1EC:
/* 003CC 80A7C1EC 24140002 */  addiu   $s4, $zero, 0x0002         ## $s4 = 00000002
/* 003D0 80A7C1F0 12540002 */  beq     $s2, $s4, .L80A7C1FC       
/* 003D4 80A7C1F4 24010003 */  addiu   $at, $zero, 0x0003         ## $at = 00000003
/* 003D8 80A7C1F8 1641002E */  bne     $s2, $at, .L80A7C2B4       
.L80A7C1FC:
/* 003DC 80A7C1FC 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 003E0 80A7C200 0C29EFE8 */  jal     func_80A7BFA0              
/* 003E4 80A7C204 02602825 */  or      $a1, $s3, $zero            ## $a1 = 00000000
/* 003E8 80A7C208 10400007 */  beq     $v0, $zero, .L80A7C228     
/* 003EC 80A7C20C 00000000 */  nop
/* 003F0 80A7C210 962D0314 */  lhu     $t5, 0x0314($s1)           ## 00000314
/* 003F4 80A7C214 44804000 */  mtc1    $zero, $f8                 ## $f8 = 0.00
/* 003F8 80A7C218 3C0180A8 */  lui     $at, %hi(D_80A7DEB0)       ## $at = 80A80000
/* 003FC 80A7C21C 35AE0010 */  ori     $t6, $t5, 0x0010           ## $t6 = 00000010
/* 00400 80A7C220 A62E0314 */  sh      $t6, 0x0314($s1)           ## 00000314
/* 00404 80A7C224 E428DEB0 */  swc1    $f8, %lo(D_80A7DEB0)($at)  
.L80A7C228:
/* 00408 80A7C228 1654001A */  bne     $s2, $s4, .L80A7C294       
/* 0040C 80A7C22C 00000000 */  nop
/* 00410 80A7C230 A6200034 */  sh      $zero, 0x0034($s1)         ## 00000034
/* 00414 80A7C234 862F0034 */  lh      $t7, 0x0034($s1)           ## 00000034
/* 00418 80A7C238 00008025 */  or      $s0, $zero, $zero          ## $s0 = 00000000
/* 0041C 80A7C23C 26721C24 */  addiu   $s2, $s3, 0x1C24           ## $s2 = 00001C24
/* 00420 80A7C240 A62F00B8 */  sh      $t7, 0x00B8($s1)           ## 000000B8
/* 00424 80A7C244 C62A0028 */  lwc1    $f10, 0x0028($s1)          ## 00000028
.L80A7C248:
/* 00428 80A7C248 8E270024 */  lw      $a3, 0x0024($s1)           ## 00000024
/* 0042C 80A7C24C 24090003 */  addiu   $t1, $zero, 0x0003         ## $t1 = 00000003
/* 00430 80A7C250 E7AA0010 */  swc1    $f10, 0x0010($sp)          
/* 00434 80A7C254 C630002C */  lwc1    $f16, 0x002C($s1)          ## 0000002C
/* 00438 80A7C258 02402025 */  or      $a0, $s2, $zero            ## $a0 = 00001C24
/* 0043C 80A7C25C 02602825 */  or      $a1, $s3, $zero            ## $a1 = 00000000
/* 00440 80A7C260 E7B00014 */  swc1    $f16, 0x0014($sp)          
/* 00444 80A7C264 863800B4 */  lh      $t8, 0x00B4($s1)           ## 000000B4
/* 00448 80A7C268 24060020 */  addiu   $a2, $zero, 0x0020         ## $a2 = 00000020
/* 0044C 80A7C26C AFB80018 */  sw      $t8, 0x0018($sp)           
/* 00450 80A7C270 863900B6 */  lh      $t9, 0x00B6($s1)           ## 000000B6
/* 00454 80A7C274 AFB9001C */  sw      $t9, 0x001C($sp)           
/* 00458 80A7C278 862800B8 */  lh      $t0, 0x00B8($s1)           ## 000000B8
/* 0045C 80A7C27C AFA90024 */  sw      $t1, 0x0024($sp)           
/* 00460 80A7C280 0C00C7D4 */  jal     Actor_Spawn
              ## ActorSpawn
/* 00464 80A7C284 AFA80020 */  sw      $t0, 0x0020($sp)           
/* 00468 80A7C288 26100001 */  addiu   $s0, $s0, 0x0001           ## $s0 = 00000001
/* 0046C 80A7C28C 5614FFEE */  bnel    $s0, $s4, .L80A7C248       
/* 00470 80A7C290 C62A0028 */  lwc1    $f10, 0x0028($s1)          ## 00000028
.L80A7C294:
/* 00474 80A7C294 0C29F4E7 */  jal     func_80A7D39C              
/* 00478 80A7C298 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 0047C 80A7C29C 3C0280A8 */  lui     $v0, %hi(D_80A7DEB8)       ## $v0 = 80A80000
/* 00480 80A7C2A0 2442DEB8 */  addiu   $v0, $v0, %lo(D_80A7DEB8)  ## $v0 = 80A7DEB8
/* 00484 80A7C2A4 844A0000 */  lh      $t2, 0x0000($v0)           ## 80A7DEB8
/* 00488 80A7C2A8 254B0001 */  addiu   $t3, $t2, 0x0001           ## $t3 = 00000001
/* 0048C 80A7C2AC 10000019 */  beq     $zero, $zero, .L80A7C314   
/* 00490 80A7C2B0 A44B0000 */  sh      $t3, 0x0000($v0)           ## 80A7DEB8
.L80A7C2B4:
/* 00494 80A7C2B4 0C03F66B */  jal     Math_Rand_ZeroOne
              ## Rand.Next() float
/* 00498 80A7C2B8 00000000 */  nop
/* 0049C 80A7C2BC 3C0180A8 */  lui     $at, %hi(D_80A7DFD4)       ## $at = 80A80000
/* 004A0 80A7C2C0 C432DFD4 */  lwc1    $f18, %lo(D_80A7DFD4)($at) 
/* 004A4 80A7C2C4 4612003C */  c.lt.s  $f0, $f18                  
/* 004A8 80A7C2C8 00000000 */  nop
/* 004AC 80A7C2CC 45000005 */  bc1f    .L80A7C2E4                 
/* 004B0 80A7C2D0 00000000 */  nop
/* 004B4 80A7C2D4 0C29F0E8 */  jal     func_80A7C3A0              
/* 004B8 80A7C2D8 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 004BC 80A7C2DC 1000000E */  beq     $zero, $zero, .L80A7C318   
/* 004C0 80A7C2E0 8FBF0044 */  lw      $ra, 0x0044($sp)           
.L80A7C2E4:
/* 004C4 80A7C2E4 3C0180A8 */  lui     $at, %hi(D_80A7DFD8)       ## $at = 80A80000
/* 004C8 80A7C2E8 C424DFD8 */  lwc1    $f4, %lo(D_80A7DFD8)($at)  
/* 004CC 80A7C2EC 4604003C */  c.lt.s  $f0, $f4                   
/* 004D0 80A7C2F0 00000000 */  nop
/* 004D4 80A7C2F4 45000005 */  bc1f    .L80A7C30C                 
/* 004D8 80A7C2F8 00000000 */  nop
/* 004DC 80A7C2FC 0C29F166 */  jal     func_80A7C598              
/* 004E0 80A7C300 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
/* 004E4 80A7C304 10000004 */  beq     $zero, $zero, .L80A7C318   
/* 004E8 80A7C308 8FBF0044 */  lw      $ra, 0x0044($sp)           
.L80A7C30C:
/* 004EC 80A7C30C 0C29F206 */  jal     func_80A7C818              
/* 004F0 80A7C310 02202025 */  or      $a0, $s1, $zero            ## $a0 = 00000000
.L80A7C314:
/* 004F4 80A7C314 8FBF0044 */  lw      $ra, 0x0044($sp)           
.L80A7C318:
/* 004F8 80A7C318 8FB00030 */  lw      $s0, 0x0030($sp)           
/* 004FC 80A7C31C 8FB10034 */  lw      $s1, 0x0034($sp)           
/* 00500 80A7C320 8FB20038 */  lw      $s2, 0x0038($sp)           
/* 00504 80A7C324 8FB3003C */  lw      $s3, 0x003C($sp)           
/* 00508 80A7C328 8FB40040 */  lw      $s4, 0x0040($sp)           
/* 0050C 80A7C32C 03E00008 */  jr      $ra                        
/* 00510 80A7C330 27BD0048 */  addiu   $sp, $sp, 0x0048           ## $sp = 00000000


