glabel EnCs_Update
/* 00B74 809E2424 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 00B78 809E2428 AFBF001C */  sw      $ra, 0x001C($sp)           
/* 00B7C 809E242C AFB00018 */  sw      $s0, 0x0018($sp)           
/* 00B80 809E2430 AFA50034 */  sw      $a1, 0x0034($sp)           
/* 00B84 809E2434 8C820210 */  lw      $v0, 0x0210($a0)           ## 00000210
/* 00B88 809E2438 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00B8C 809E243C 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 00B90 809E2440 1440000F */  bne     $v0, $zero, .L809E2480     
/* 00B94 809E2444 00000000 */  nop
/* 00B98 809E2448 C4840164 */  lwc1    $f4, 0x0164($a0)           ## 00000164
/* 00B9C 809E244C 24010009 */  addiu   $at, $zero, 0x0009         ## $at = 00000009
/* 00BA0 809E2450 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00BA4 809E2454 4600218D */  trunc.w.s $f6, $f4                   
/* 00BA8 809E2458 44023000 */  mfc1    $v0, $f6                   
/* 00BAC 809E245C 00000000 */  nop
/* 00BB0 809E2460 10410003 */  beq     $v0, $at, .L809E2470       
/* 00BB4 809E2464 24010017 */  addiu   $at, $zero, 0x0017         ## $at = 00000017
/* 00BB8 809E2468 54410023 */  bnel    $v0, $at, .L809E24F8       
/* 00BBC 809E246C 26060194 */  addiu   $a2, $s0, 0x0194           ## $a2 = 00000194
.L809E2470:
/* 00BC0 809E2470 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 00BC4 809E2474 240528EF */  addiu   $a1, $zero, 0x28EF         ## $a1 = 000028EF
/* 00BC8 809E2478 1000001F */  beq     $zero, $zero, .L809E24F8   
/* 00BCC 809E247C 26060194 */  addiu   $a2, $s0, 0x0194           ## $a2 = 00000194
.L809E2480:
/* 00BD0 809E2480 54410010 */  bnel    $v0, $at, .L809E24C4       
/* 00BD4 809E2484 24010002 */  addiu   $at, $zero, 0x0002         ## $at = 00000002
/* 00BD8 809E2488 C6080164 */  lwc1    $f8, 0x0164($s0)           ## 00000164
/* 00BDC 809E248C 2401000A */  addiu   $at, $zero, 0x000A         ## $at = 0000000A
/* 00BE0 809E2490 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00BE4 809E2494 4600428D */  trunc.w.s $f10, $f8                  
/* 00BE8 809E2498 44025000 */  mfc1    $v0, $f10                  
/* 00BEC 809E249C 00000000 */  nop
/* 00BF0 809E24A0 10410003 */  beq     $v0, $at, .L809E24B0       
/* 00BF4 809E24A4 24010019 */  addiu   $at, $zero, 0x0019         ## $at = 00000019
/* 00BF8 809E24A8 54410013 */  bnel    $v0, $at, .L809E24F8       
/* 00BFC 809E24AC 26060194 */  addiu   $a2, $s0, 0x0194           ## $a2 = 00000194
.L809E24B0:
/* 00C00 809E24B0 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 00C04 809E24B4 240528EF */  addiu   $a1, $zero, 0x28EF         ## $a1 = 000028EF
/* 00C08 809E24B8 1000000F */  beq     $zero, $zero, .L809E24F8   
/* 00C0C 809E24BC 26060194 */  addiu   $a2, $s0, 0x0194           ## $a2 = 00000194
/* 00C10 809E24C0 24010002 */  addiu   $at, $zero, 0x0002         ## $at = 00000002
.L809E24C4:
/* 00C14 809E24C4 5441000C */  bnel    $v0, $at, .L809E24F8       
/* 00C18 809E24C8 26060194 */  addiu   $a2, $s0, 0x0194           ## $a2 = 00000194
/* 00C1C 809E24CC C6100164 */  lwc1    $f16, 0x0164($s0)          ## 00000164
/* 00C20 809E24D0 24010014 */  addiu   $at, $zero, 0x0014         ## $at = 00000014
/* 00C24 809E24D4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00C28 809E24D8 4600848D */  trunc.w.s $f18, $f16                 
/* 00C2C 809E24DC 44199000 */  mfc1    $t9, $f18                  
/* 00C30 809E24E0 00000000 */  nop
/* 00C34 809E24E4 57210004 */  bnel    $t9, $at, .L809E24F8       
/* 00C38 809E24E8 26060194 */  addiu   $a2, $s0, 0x0194           ## $a2 = 00000194
/* 00C3C 809E24EC 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 00C40 809E24F0 240528EF */  addiu   $a1, $zero, 0x28EF         ## $a1 = 000028EF
/* 00C44 809E24F4 26060194 */  addiu   $a2, $s0, 0x0194           ## $a2 = 00000194
.L809E24F8:
/* 00C48 809E24F8 00C02825 */  or      $a1, $a2, $zero            ## $a1 = 00000194
/* 00C4C 809E24FC AFA60024 */  sw      $a2, 0x0024($sp)           
/* 00C50 809E2500 0C0189B7 */  jal     ActorCollider_Cylinder_Update
              
/* 00C54 809E2504 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00C58 809E2508 8FA40034 */  lw      $a0, 0x0034($sp)           
/* 00C5C 809E250C 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00C60 809E2510 34211E60 */  ori     $at, $at, 0x1E60           ## $at = 00011E60
/* 00C64 809E2514 8FA60024 */  lw      $a2, 0x0024($sp)           
/* 00C68 809E2518 0C017713 */  jal     Actor_CollisionCheck_SetOT
              ## CollisionCheck_setOT
/* 00C6C 809E251C 00812821 */  addu    $a1, $a0, $at              
/* 00C70 809E2520 8E190190 */  lw      $t9, 0x0190($s0)           ## 00000190
/* 00C74 809E2524 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00C78 809E2528 8FA50034 */  lw      $a1, 0x0034($sp)           
/* 00C7C 809E252C 0320F809 */  jalr    $ra, $t9                   
/* 00C80 809E2530 00000000 */  nop
/* 00C84 809E2534 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00C88 809E2538 0C27874E */  jal     func_809E1D38              
/* 00C8C 809E253C 8FA50034 */  lw      $a1, 0x0034($sp)           
/* 00C90 809E2540 8E0901E8 */  lw      $t1, 0x01E8($s0)           ## 000001E8
/* 00C94 809E2544 3C0F809E */  lui     $t7, %hi(D_809E2958)       ## $t7 = 809E0000
/* 00C98 809E2548 252AFFFF */  addiu   $t2, $t1, 0xFFFF           ## $t2 = FFFFFFFF
/* 00C9C 809E254C 0541000C */  bgez    $t2, .L809E2580            
/* 00CA0 809E2550 AE0A01E8 */  sw      $t2, 0x01E8($s0)           ## 000001E8
/* 00CA4 809E2554 8E0C01E4 */  lw      $t4, 0x01E4($s0)           ## 000001E4
/* 00CA8 809E2558 25820001 */  addiu   $v0, $t4, 0x0001           ## $v0 = 00000001
/* 00CAC 809E255C 28410003 */  slti    $at, $v0, 0x0003           
/* 00CB0 809E2560 14200003 */  bne     $at, $zero, .L809E2570     
/* 00CB4 809E2564 AE0201E4 */  sw      $v0, 0x01E4($s0)           ## 000001E4
/* 00CB8 809E2568 AE0001E4 */  sw      $zero, 0x01E4($s0)         ## 000001E4
/* 00CBC 809E256C 00001025 */  or      $v0, $zero, $zero          ## $v0 = 00000000
.L809E2570:
/* 00CC0 809E2570 00027080 */  sll     $t6, $v0,  2               
/* 00CC4 809E2574 01EE7821 */  addu    $t7, $t7, $t6              
/* 00CC8 809E2578 8DEF2958 */  lw      $t7, %lo(D_809E2958)($t7)  
/* 00CCC 809E257C AE0F01E8 */  sw      $t7, 0x01E8($s0)           ## 000001E8
.L809E2580:
/* 00CD0 809E2580 8FBF001C */  lw      $ra, 0x001C($sp)           
/* 00CD4 809E2584 8FB00018 */  lw      $s0, 0x0018($sp)           
/* 00CD8 809E2588 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000
/* 00CDC 809E258C 03E00008 */  jr      $ra                        
/* 00CE0 809E2590 00000000 */  nop
