glabel func_80B86494
/* 00364 80B86494 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 00368 80B86498 3C0180B8 */  lui     $at, %hi(D_80B8683C)       ## $at = 80B80000
/* 0036C 80B8649C C424683C */  lwc1    $f4, %lo(D_80B8683C)($at)  
/* 00370 80B864A0 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 00374 80B864A4 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 00378 80B864A8 3C01C0A0 */  lui     $at, 0xC0A0                ## $at = C0A00000
/* 0037C 80B864AC 44813000 */  mtc1    $at, $f6                   ## $f6 = -5.00
/* 00380 80B864B0 44800000 */  mtc1    $zero, $f0                 ## $f0 = 0.00
/* 00384 80B864B4 3C0140C0 */  lui     $at, 0x40C0                ## $at = 40C00000
/* 00388 80B864B8 44814000 */  mtc1    $at, $f8                   ## $f8 = 6.00
/* 0038C 80B864BC 3C0580B8 */  lui     $a1, %hi(func_80B862EC)    ## $a1 = 80B80000
/* 00390 80B864C0 24A562EC */  addiu   $a1, $a1, %lo(func_80B862EC) ## $a1 = 80B862EC
/* 00394 80B864C4 E484006C */  swc1    $f4, 0x006C($a0)           ## 0000006C
/* 00398 80B864C8 E4860070 */  swc1    $f6, 0x0070($a0)           ## 00000070
/* 0039C 80B864CC E480005C */  swc1    $f0, 0x005C($a0)           ## 0000005C
/* 003A0 80B864D0 E4800064 */  swc1    $f0, 0x0064($a0)           ## 00000064
/* 003A4 80B864D4 0C2E184C */  jal     func_80B86130              
/* 003A8 80B864D8 E4880060 */  swc1    $f8, 0x0060($a0)           ## 00000060
/* 003AC 80B864DC 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 003B0 80B864E0 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 003B4 80B864E4 03E00008 */  jr      $ra                        
/* 003B8 80B864E8 00000000 */  nop


