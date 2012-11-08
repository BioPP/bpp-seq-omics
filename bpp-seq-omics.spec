%define _basename bpp-seq-omics
%define _version 2.0.3
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
BuildRequires: cmake >= 2.6.0
BuildRequires: gcc-c++ >= 4.0.0
BuildRequires: libbpp-core2 = %{_version}
BuildRequires: libbpp-core-devel = %{_version}
BuildRequires: libbpp-seq9 = %{_version}
BuildRequires: libbpp-seq-devel = %{_version}

AutoReq: yes
AutoProv: yes

%description
This library contains the genomics components of the Bio++ sequence library.
It is part of the Bio++ project.

%package -n libbpp-seq-omics1
Summary: Bio++ Sequence library: genomics components
Group: Development/Libraries/C and C++

%description -n libbpp-seq-omics1
This library contains the genomics components of the Bio++ sequence library.
It is part of the Bio++ project.


%package -n libbpp-seq-omics-devel
Summary: Bio++ Sequence library: genomics components
Group: Development/Libraries/C and C++
Requires: libbpp-seq-omics1 = %{_version}
Requires: libbpp-seq9 = %{_version}
Requires: libbpp-seq-devel = %{_version}
Requires: libbpp-core2 = %{_version}
Requires: libbpp-core-devel = %{_version}

%description -n libbpp-seq-omics-devel
The libbpp-seq-omics-devel package contains the header files and static libraries for
building applications which use %{_basename}.

%prep
%setup -q

%build
CFLAGS="$RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix} -DBUILD_TESTING=OFF"
if [ %{_lib} == 'lib64' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DLIB_SUFFIX=64"
fi
cmake $CMAKE_FLAGS .
make

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -n libbpp-seq-omics1 -p /sbin/ldconfig

%postun -n libbpp-seq-omics1 -p /sbin/ldconfig

%files -n libbpp-seq-omics1
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/%{_lib}/lib*.so.*

%files -n libbpp-seq-omics-devel
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/%{_lib}/lib*.so
%{_prefix}/%{_lib}/lib*.a
%{_prefix}/include/*

%changelog
* Tue Nov 06 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.0.3-1
- First draft of the spec file

