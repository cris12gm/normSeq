from django.contrib import admin
from django.urls import include, re_path
from django.conf import settings
from django.conf.urls.static import static
from .views import statusJob

urlpatterns = [
    re_path('', statusJob.as_view(), name='status'),
    ]+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)