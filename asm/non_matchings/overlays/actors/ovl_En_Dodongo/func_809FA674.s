glabel func_809FA674
/* 02424 809FA674 27BDFFE8 */  addiu   $sp, $sp, 0xFFE8           ## $sp = FFFFFFE8
/* 02428 809FA678 AFBF0014 */  sw      $ra, 0x0014($sp)           
/* 0242C 809FA67C AFA40018 */  sw      $a0, 0x0018($sp)           
/* 02430 809FA680 AFA5001C */  sw      $a1, 0x001C($sp)           
/* 02434 809FA684 AFA60020 */  sw      $a2, 0x0020($sp)           
/* 02438 809FA688 0C01DE1C */  jal     Math_Sins
              ## sins?
/* 0243C 809FA68C 87A4001A */  lh      $a0, 0x001A($sp)           
/* 02440 809FA690 C7A6001C */  lwc1    $f6, 0x001C($sp)           
/* 02444 809FA694 8FA20020 */  lw      $v0, 0x0020($sp)           
/* 02448 809FA698 46060202 */  mul.s   $f8, $f0, $f6              
/* 0244C 809FA69C C4440000 */  lwc1    $f4, 0x0000($v0)           ## 00000000
/* 02450 809FA6A0 46082280 */  add.s   $f10, $f4, $f8             
/* 02454 809FA6A4 E44A0000 */  swc1    $f10, 0x0000($v0)          ## 00000000
/* 02458 809FA6A8 0C01DE0D */  jal     Math_Coss
              ## coss?
/* 0245C 809FA6AC 87A4001A */  lh      $a0, 0x001A($sp)           
/* 02460 809FA6B0 C7B2001C */  lwc1    $f18, 0x001C($sp)          
/* 02464 809FA6B4 8FA20020 */  lw      $v0, 0x0020($sp)           
/* 02468 809FA6B8 46120182 */  mul.s   $f6, $f0, $f18             
/* 0246C 809FA6BC C4500008 */  lwc1    $f16, 0x0008($v0)          ## 00000008
/* 02470 809FA6C0 46068100 */  add.s   $f4, $f16, $f6             
/* 02474 809FA6C4 E4440008 */  swc1    $f4, 0x0008($v0)           ## 00000008
/* 02478 809FA6C8 8FBF0014 */  lw      $ra, 0x0014($sp)           
/* 0247C 809FA6CC 27BD0018 */  addiu   $sp, $sp, 0x0018           ## $sp = 00000000
/* 02480 809FA6D0 03E00008 */  jr      $ra                        
/* 02484 809FA6D4 00000000 */  nop


