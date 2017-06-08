%define _basename bpp-seq-omics
%define _version 2.3.1
%define _release 1
%define _prefix /usr

URL: http://biopp.univ-montp2.fr/

Name: %{_basename}
Version: %{_version}
Release: %{_release}
License: CECILL-2.0
Vendor: The Bio++ Project
Source: http://biopp.univ-montp2.fr/repos/sources/%{_basename}-%{_version}.tar.gz
Summary: Bio++ Sequence library: genomics components
Group: Development/Libraries/C and C++
Requires: bpp-core = %{_version}
Requires: bpp-seq = %{_version}

BuildRoot: %{_builddir}/%{_basename}-root
BuildRequires: cmake >= 2.8.11
BuildRequires: gcc-c++ >= 4.7.0
BuildRequires: libbpp-core3 = %{_version}
BuildRequires: libbpp-core-devel = %{_version}
BuildRequires: libbpp-seq11 = %{_version}
BuildRequires: libbpp-seq-devel = %{_version}

AutoReq: yes
AutoProv: yes

%description
This library contains the genomics components of the Bio++ sequence library.
It is part of the Bio++ project.

%package -n libbpp-seq-omics2
Summary: Bio++ Sequence library: genomics components
Group: Development/Libraries/C and C++

%description -n libbpp-seq-omics2
This library contains the genomics components of the Bio++ sequence library.
It is part of the Bio++ project.

%package -n libbpp-seq-omics-devel
Summary: Bio++ Sequence library: genomics components
Group: Development/Libraries/C and C++
Requires: libbpp-seq-omics2 = %{_version}
Requires: libbpp-seq11 = %{_version}
Requires: libbpp-seq-devel = %{_version}
Requires: libbpp-core3 = %{_version}
Requires: libbpp-core-devel = %{_version}

%description -n libbpp-seq-omics-devel
The libbpp-seq-omics-devel package contains the header files and static libraries for
building applications which use %{_basename}.

%prep
%setup -q

%build
CFLAGS="$RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix} -DBUILD_TESTING=OFF"
cmake $CMAKE_FLAGS .
make

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -n libbpp-seq-omics2 -p /sbin/ldconfig

%postun -n libbpp-seq-omics2 -p /sbin/ldconfig

%files -n libbpp-seq-omics2
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/%{_lib}/lib*.so.*

%files -n libbpp-seq-omics-devel
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%dir %{_prefix}/%{_lib}/cmake/
%dir %{_prefix}/%{_lib}/cmake/bpp-seq-omics
%{_prefix}/%{_lib}/lib*.so
%{_prefix}/%{_lib}/lib*.a
%{_prefix}/%{_lib}/cmake/bpp-seq-omics/bpp-seq-omics*.cmake
%{_prefix}/include/*

%changelog
* Tue Jun 06 2017 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.3.1-1
- Increased interface number
* Wed May 10 2017 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.3.0-1
- Several bugs fixed and performance improvements
- New Maf filters: LiftOver, RemoveEmptySequence
- Support for BED features input
- New output options: MSMC, PLINK
- More options for Vcf and Alignment output
- Upgrade to C++11
* Mon Sep 22 2014 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.2.0-1
- New statistics, including sequence diversity estimators
- Several bugs and memory leaks fixed.
* Wed Mar 06 2013 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.1.0-1
- Maf to VCF tool added as a MafIterator.
* Tue Nov 06 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.0.3-1
- First draft of the spec file

