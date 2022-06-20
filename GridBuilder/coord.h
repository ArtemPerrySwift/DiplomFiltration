#pragma once
#define USER_COORDS true

struct Coord
{
	double x;
	double y;
	double z;

	Coord() { x = y = z = 0; }

	Coord(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	void transIntoMas(double mas[3])
	{
		mas[0] = x;
		mas[1] = y;
		mas[2] = z;
	}

	Coord operator +(const Coord& coord2)
	{
		Coord newCoord;
		newCoord.x = this->x + coord2.x;
		newCoord.y = this->y + coord2.y;
		newCoord.z = this->z + coord2.z;
		return newCoord;
	}

	Coord operator -(const Coord& coord2)
	{
		Coord newCoord;
		newCoord.x = this->x - coord2.x;
		newCoord.y = this->y - coord2.y;
		newCoord.z = this->z - coord2.z;
		return newCoord;
	}

	Coord operator /(const Coord& coord2)
	{
		Coord newCoord;
		newCoord.x = this->x / coord2.x;
		newCoord.y = this->y / coord2.y;
		newCoord.z = this->z / coord2.z;
		return newCoord;
	}

	Coord operator *(const Coord& coord2)
	{
		Coord newCoord;
		newCoord.x = this->x * coord2.x;
		newCoord.y = this->y * coord2.y;
		newCoord.z = this->z * coord2.z;
		return newCoord;
	}

	Coord operator /(const double& a)
	{
		Coord newCoord;
		newCoord.x = this->x / a;
		newCoord.y = this->y / a;
		newCoord.z = this->z / a;
		return newCoord;
	}

	Coord operator *(const double& a)
	{
		Coord newCoord;
		newCoord.x = this->x * a;
		newCoord.y = this->y * a;
		newCoord.z = this->z * a;
		return newCoord;
	}

	Coord& operator +=(const Coord& coord2)
	{
		Coord newCoord;
		this->x += coord2.x;
		this->y += coord2.y;
		this->z += coord2.z;
		return *this;
	}

	Coord& operator *=(const Coord& coord2)
	{
		Coord newCoord;
		this->x *= coord2.x;
		this->y *= coord2.y;
		this->z *= coord2.z;
		return *this;
	}

	Coord& operator *=(const double& a)
	{
		Coord newCoord;
		this->x *= a;
		this->y *= a;
		this->z *= a;
		return *this;
	}

	Coord& operator =(const double a[3])
	{
		Coord newCoord;
		this->x = a[0];
		this->y = a[1];
		this->z = a[2];
		return *this;
	}
};

namespace ecoord
{
	struct ECoord
	{
		double ksi;
		double nu;
		double etta;

		ECoord() { ksi = nu = etta = 0; }

		ECoord(double x, double y, double z)
		{
			this->ksi = x;
			this->nu = y;
			this->etta = z;
		}

		void transIntoMas(double mas[3])
		{
			mas[0] = ksi;
			mas[1] = nu;
			mas[2] = etta;
		}

		ECoord operator +(const ECoord& coord2)
		{
			ECoord newCoord;
			newCoord.ksi = this->ksi + coord2.ksi;
			newCoord.nu = this->nu + coord2.nu;
			newCoord.etta = this->etta + coord2.etta;
			return newCoord;
		}

		ECoord operator -(const ECoord& coord2)
		{
			ECoord newCoord;
			newCoord.ksi = this->ksi - coord2.ksi;
			newCoord.nu = this->nu - coord2.nu;
			newCoord.etta = this->etta - coord2.etta;
			return newCoord;
		}

		ECoord operator /(const ECoord& coord2)
		{
			ECoord newCoord;
			newCoord.ksi = this->ksi / coord2.ksi;
			newCoord.nu = this->nu / coord2.nu;
			newCoord.etta = this->etta / coord2.etta;
			return newCoord;
		}

		ECoord operator *(const ECoord& coord2)
		{
			ECoord newCoord;
			newCoord.ksi = this->ksi * coord2.ksi;
			newCoord.nu = this->nu * coord2.nu;
			newCoord.etta = this->etta * coord2.etta;
			return newCoord;
		}

		ECoord operator /(const double& a)
		{
			ECoord newCoord;
			newCoord.ksi = this->ksi / a;
			newCoord.nu = this->nu / a;
			newCoord.etta = this->etta / a;
			return newCoord;
		}

		ECoord operator *(const double& a)
		{
			ECoord newCoord;
			newCoord.ksi = this->ksi * a;
			newCoord.nu = this->nu * a;
			newCoord.etta = this->etta * a;
			return newCoord;
		}

		ECoord& operator +=(const ECoord& coord2)
		{
			ECoord newCoord;
			this->ksi += coord2.ksi;
			this->nu += coord2.nu;
			this->etta += coord2.etta;
			return *this;
		}

		ECoord& operator *=(const ECoord& coord2)
		{
			ECoord newCoord;
			this->ksi *= coord2.ksi;
			this->nu *= coord2.nu;
			this->etta *= coord2.etta;
			return *this;
		}

		ECoord& operator *=(const double& a)
		{
			ECoord newCoord;
			this->ksi *= a;
			this->nu *= a;
			this->etta *= a;
			return *this;
		}

		ECoord& operator =(const double a[3])
		{
			ECoord newCoord;
			this->ksi = a[0];
			this->nu = a[1];
			this->etta = a[2];
			return *this;
		}

		ECoord& operator =(const Coord& p)
		{
			this->ksi = p.x;
			this->nu = p.y;
			this->etta = p.z;
			return *this;
		}

		operator Coord()
		{
			Coord p(this->ksi, this->nu, this->etta);
			return p;
		}
	};
	
}
