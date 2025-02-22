{
  description = "ISPACK3";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
   { nixpkgs, flake-utils, ... }:
   flake-utils.lib.eachDefaultSystem (
     system:
     let
       pkgs = nixpkgs.legacyPackages.${system};
       ispack3 = pkgs.stdenv.mkDerivation {
         pname = "ispack3";
         version = "3.2.2";
         src = ./ispack-3.2.2;
         nativeBuildInputs = with pkgs; [ gfortran mpi gnumake ];
         buildPhase = ''
           make
           make clean
         '';
         installPhase = ''
           mkdir -p $out/lib
           cp ./libispack3.a $out/lib/libispack3.a
         '';
       };
     in
     {
       packages = {
         inherit ispack3;
         default = ispack3;
       };
     }
   );
}
